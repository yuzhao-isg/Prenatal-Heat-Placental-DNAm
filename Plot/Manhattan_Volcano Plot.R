# Manhattan and volcano plots for all participants, females, and males.
rm(list = ls())
library(ggplot2)
library(dplyr)
library(ggrepel)
library(ggnewscale)

# Paths and settings ------------------------------------------------------
results_dir <- file.path("~/010_EWAS_temp_YZ/results") # # change to your path
qc_dir <- file.path(results_dir, "Main_Model/CpG_quality_control")
annotation_dir <- file.path(results_dir, "Downstream/Annotation")
plot_dir <- file.path(results_dir, "Plot")

suggestive_p <- 1e-5
bonferroni_p <- 0.05 / 688634
suggestive_y <- -log10(suggestive_p)
bonferroni_y <- -log10(bonferroni_p)
chromosomes <- c(as.character(1:22), "X", "Y")
exposure_colors <- c(
  "Indoor Temperature" = "#8B0000", "Indoor Heat Index" = "#F08080",
  "Outdoor Temperature" = "#6A0DAD", "Outdoor Heat Index" = "#C49FFF"
)
exposures <- data.frame(
  location = c("Indoor", "Indoor", "Outdoor", "Outdoor"),
  code = c("indoor_HI", "indoor_T", "outdoor_HI", "outdoor_T"),
  metric = c("HI", "heat", "HI", "heat"),
  period = c("indoor_day", "indoor_day", "home", "home"),
  Exposure = c("Indoor Heat Index", "Indoor Temperature",
               "Outdoor Heat Index", "Outdoor Temperature"),
  stringsAsFactors = FALSE
)

# Input and gene-name helpers ---------------------------------------------
read_rdata <- function(path, object) {
  # Windows' extended UNC prefix allows long network paths.
  if (.Platform$OS.type == "windows" && startsWith(path, "//")) {
    path <- gsub("/+", "/", substring(path, 3))
    path <- paste0("\\\\?\\UNC\\", chartr("/", "\\", path))
  }
  data <- new.env()
  load(path, envir = data)
  data[[object]]
}

unique_gene_names <- function(values) {
  vapply(as.character(values), function(value) {
    if (is.na(value)) return(NA_character_)
    genes <- trimws(strsplit(value, ";", fixed = TRUE)[[1]])
    paste(unique(genes[nzchar(genes)]), collapse = ";")
  }, character(1), USE.NAMES = FALSE)
}

read_group <- function(sex = "") {
  suffix <- if (nzchar(sex)) paste0("_", sex) else ""
  rows <- lapply(seq_len(nrow(exposures)), function(i) {
    exposure <- exposures[i, ]
    directory <- file.path(qc_dir, exposure$location)
    if (nzchar(sex)) {
      directory <- file.path(directory, sprintf(
        "count_%s_exposure_95_%s", exposure$metric, exposure$period
      ))
    }
    filename <- paste0("Exclude_CpGs_results_", exposure$code, "95", suffix, ".RData")
    object <- paste0("result_", exposure$code, "95", suffix)
    data <- read_rdata(file.path(directory, filename), object)
    exposure_name <- if (nzchar(sex)) paste(exposure$Exposure, sex) else exposure$Exposure
    mutate(data, Exposure = exposure_name)
  })
  do.call(rbind, rows)
}

read_annotation <- function(location) {
  read.csv(file.path(annotation_dir, location, "Sig_DMPs_annotated.csv")) %>%
    rename(Suggestive = Suggested, Suggestive_4 = Suggested_4)
}

# Prepare annotations once for all three groups.
ewas_sig <- rbind(read_annotation("Indoor"), read_annotation("Outdoor")) %>%
  filter(!grepl("_M$|_F$", Exposure))
# Retain the same annotation columns as the original script.
ewas_sig <- ewas_sig[, -c(1, 3:33)]

prepare_plot_data <- function(CpGs_all) {
  data <- left_join(CpGs_all, ewas_sig, by = c("CpG", "Exposure")) %>%
    mutate(
      UCSC_RefGene_Name = unique_gene_names(UCSC_RefGene_Name),
      label = ifelse(
        Bonferroni == "yes",
        ifelse(!is.na(UCSC_RefGene_Name) & UCSC_RefGene_Name != "",
               UCSC_RefGene_Name, CpG),
        NA
      )
    )

  data$CpG_chrm <- as.character(data$CpG_chrm)
  data <- data %>%
    mutate(chr = factor(gsub("chr", "", CpG_chrm), levels = chromosomes),
           pos = (CpG_beg + CpG_end) / 2) %>%
    arrange(chr, pos)
  chr_lengths <- data %>%
    group_by(chr) %>%
    summarise(max_bp = max(pos), .groups = "drop") %>%
    arrange(chr) %>%
    mutate(cumlen = cumsum(as.numeric(max_bp)) - as.numeric(max_bp))
  chr_lengths$chr <- as.character(chr_lengths$chr)
  data <- data %>%
    left_join(select(chr_lengths, chr, cumlen), by = "chr") %>%
    mutate(BPcum = pos + cumlen)

  data %>%
    mutate(
      Exposure = trimws(as.character(Exposure)),
      significance = case_when(
        pval < bonferroni_p ~ "Bonferroni",
        pval < suggestive_p ~ "Suggestive",
        TRUE ~ NA_character_
      ),
      significance = factor(significance, levels = c("Suggestive", "Bonferroni")),
      chr_index = match(gsub("^chr", "", CpG_chrm), chromosomes),
      chr_color = if_else(chr_index %% 2 == 0, "grey20", "grey70"),
      plot_label = case_when(
        significance == "Bonferroni" & !is.na(UCSC_RefGene_Name) &
          UCSC_RefGene_Name != "" ~ UCSC_RefGene_Name,
        significance == "Bonferroni" ~ CpG,
        TRUE ~ NA_character_
      )
    )
}

# Shared axes, legends, and theme ------------------------------------------
plot_style <- function(x_label) {
  list(
    scale_y_continuous(
      breaks = c(0, 2, 4, suggestive_y, 6, bonferroni_y, 8),
      labels = c("0", "2", "4", "Suggestive", "6", "Bonferroni", "8")
    ),
    labs(x = x_label, y = expression(-log[10](italic(P))),
         color = "Exposure", shape = "CpG significance"),
    guides(
      color = guide_legend(order = 1, override.aes = list(shape = 16, size = 3)),
      shape = guide_legend(order = 2, override.aes = list(color = "black", size = 3))
    ),
    theme_bw(),
    theme(panel.grid = element_blank(), legend.position = "right")
  )
}

make_manhattan <- function(data, colors) {
  axis_set <- data %>%
    group_by(chr) %>%
    summarize(center = (max(BPcum) + min(BPcum)) / 2)
  ggplot(data, aes(BPcum, -log10(pval))) +
    geom_point(aes(color = chr_color), size = 1.5, alpha = 0.8, na.rm = TRUE) +
    scale_color_identity(guide = "none") +
    ggnewscale::new_scale_color() +
    geom_point(
      data = filter(data, !is.na(significance)),
      aes(color = Exposure, shape = significance),
      size = 2.8, alpha = 1, na.rm = TRUE
    ) +
    geom_text_repel(
      data = filter(data, significance == "Bonferroni", !is.na(plot_label)),
      aes(label = plot_label, color = Exposure),
      size = 3, fontface = "italic", max.overlaps = 30,
      show.legend = FALSE, seed = 123, na.rm = TRUE
    ) +
    geom_hline(yintercept = c(suggestive_y, bonferroni_y),
               color = "black", linetype = "dashed") +
    scale_x_continuous(breaks = axis_set$center, labels = axis_set$chr) +
    scale_color_manual(values = colors, breaks = names(colors)) +
    scale_shape_manual(
      values = c("Suggestive" = 15, "Bonferroni" = 17),
      breaks = c("Suggestive", "Bonferroni"), drop = FALSE
    ) +
    plot_style("Chromosome")
}

make_volcano <- function(data, colors) {
  ggplot(data, aes(x = coef * 10, y = -log10(pval))) +
    geom_point(color = "grey50", alpha = 0.6, size = 1.2, na.rm = TRUE) +
    geom_point(
      data = filter(data, !is.na(significance)),
      aes(color = Exposure, shape = significance),
      size = 2.5, alpha = 1, na.rm = TRUE
    ) +
    geom_text_repel(
      data = filter(data, significance == "Bonferroni", !is.na(label), label != ""),
      aes(label = label, color = Exposure),
      size = 3, max.overlaps = 40, show.legend = FALSE,
      fontface = "italic", seed = 123, na.rm = TRUE
    ) +
    geom_hline(yintercept = c(suggestive_y, bonferroni_y),
               linetype = "dashed", color = "black") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    scale_color_manual(values = colors, limits = names(colors), drop = FALSE) +
    scale_shape_manual(values = c("Suggestive" = 15, "Bonferroni" = 17), drop = FALSE) +
    plot_style("Effect estimate per 10-unit increase")
}

# Generate the same six output files ---------------------------------------
for (sex in c("", "F", "M")) {
  suffix <- if (nzchar(sex)) paste0("_", sex) else ""
  colors <- exposure_colors
  if (nzchar(sex)) names(colors) <- paste(names(colors), sex)

  ewas_plot <- prepare_plot_data(read_group(sex))
  print(count(ewas_plot, significance, Exposure))
  p_manhattan <- make_manhattan(ewas_plot, colors)
  p_volcano <- make_volcano(ewas_plot, colors)
  ggsave(file.path(plot_dir, paste0("manhattan_plot", suffix, ".png")),
         plot = p_manhattan, width = 12, height = 6, dpi = 600)
  ggsave(file.path(plot_dir, paste0("volcano_plot", suffix, ".png")),
         plot = p_volcano, width = 12, height = 6, dpi = 600)
}
