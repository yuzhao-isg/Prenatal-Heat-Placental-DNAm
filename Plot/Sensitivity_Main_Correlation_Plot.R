# Compare effect sizes: main vs sensitivity models, and males vs females.
rm(list = ls())
library(ggplot2)
library(dplyr)

# Paths -------------------------------------------------------------------
results_dir <- file.path("~/010_EWAS_temp_YZ/results") # change to your path
qc_dir <- file.path(results_dir, "Main_Model/CpG_quality_control")
output_dir <- file.path(results_dir, "Sensitivity_Model/CpG_quality_control")
plot_dir <- file.path(output_dir, "Plot")
annotation_dir <- file.path(results_dir, "Heat_days/Downstream/Annotation")

exposures <- data.frame(
  location = c("Indoor", "Indoor", "Outdoor", "Outdoor"),
  code = c("indoor_HI", "indoor_T", "outdoor_HI", "outdoor_T"),
  metric = c("HI", "heat", "HI", "heat"),
  period = c("indoor_day", "indoor_day", "home", "home"),
  Exposure = c("Indoor Heat Index", "Indoor Temperature",
               "Outdoor Heat Index", "Outdoor Temperature"),
  stringsAsFactors = FALSE
)

# Input and preparation helpers -------------------------------------------
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

read_qc <- function(i, threshold = 95, sex = NULL) {
  exposure <- exposures[i, ]
  object <- paste0("result_", exposure$code, threshold)
  directory <- file.path(qc_dir, exposure$location)
  if (!is.null(sex)) {
    object <- paste0(object, "_", sex)
    directory <- file.path(directory, sprintf(
      "count_%s_exposure_95_%s", exposure$metric, exposure$period
    ))
  }
  filename <- paste0("Exclude_CpGs_", sub("^result_", "results_", object), ".RData")
  read_rdata(file.path(directory, filename), object)
}

read_sensitivity <- function(i, analysis) {
  exposure <- exposures[i, ]
  stage <- if (analysis == "celltype") "Main" else "Sensitivity"
  date <- if (analysis == "celltype") "20251223" else "20260406"
  tag <- switch(analysis, ap = "AP", celltype = "MAIN", analysis)
  metric <- exposure$metric
  # Preserve the original indoor-temperature inputs for AP and GWAS ethnicity.
  if (i == 2 && analysis %in% c("ap", "gwase")) metric <- "HI"

  stem <- sprintf("count_%s_exposure_95_%s_%s", metric, exposure$period, tag)
  directory <- file.path(paste0(exposure$location, "_", stage, "_Whole_pregnancy"),
    paste0("BISC_", date, "_Output"), stem
  )
  if (analysis == "EUR") directory <- file.path(directory, "ethnicity_c_2cat_eur")
  filename <- paste0("BISC_", date, "_", stem, "_allanalyses.RData")
  model <- if (analysis == "celltype") "Adjusted" else "AdjustedwithCellType"
  read_rdata(file.path(directory, filename), "alldataout")[[model]]
}

rename_effects <- function(data, suffix) {
  columns <- c("coef", "se", "pval")
  names(data)[match(columns, names(data))] <- paste0(columns, "_", suffix)
  data
}

prepare_effects <- function(data, exposure, suffix, main = FALSE) {
  data <- mutate(data, Exposure = exposure)
  if (main) {
    data <- data %>%
      subset(Suggestive == "yes") %>%
      mutate(label = ifelse(Bonferroni == "yes", "Bonferroni", "Suggestive"))
  }
  columns <- c("CpG", "coef", "se", "pval", if (main) "label", "Exposure")
  data <- data %>% select(all_of(columns)) %>% mutate(coef = coef * 10)
  rename_effects(data, suffix)
}

read_effects <- function(analysis) {
  rows <- lapply(seq_len(nrow(exposures)), function(i) {
    if (analysis %in% c("main", "90th", "99th")) {
      threshold <- switch(analysis, main = 95, "90th" = 90, "99th" = 99)
      data <- read_qc(i, threshold)
    } else {
      data <- read_sensitivity(i, analysis)
    }
    prepare_effects(data, exposures$Exposure[i], analysis, main = analysis == "main")
  })
  do.call(rbind, rows)
}

# Shared plot styling ------------------------------------------------------
comparison_style <- function() {
  list(
    scale_x_continuous(breaks = seq(-0.1, 0.05, 0.025)),
    scale_y_continuous(breaks = seq(-0.1, 0.05, 0.025)),
    scale_color_manual(
      values = c("Indoor Temperature" = "#8B0000", "Indoor Heat Index" = "#F08080",
                 "Outdoor Temperature" = "#6A0DAD", "Outdoor Heat Index" = "#C49FFF"),
      drop = FALSE
    ),
    theme_classic(),
    theme(
      axis.line = element_line(color = "white"),
      axis.ticks = element_line(color = "white"),
      panel.grid = element_blank(),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10)
    )
  )
}

plot_sensitivity <- function(data, analysis) {
  ggplot(data, aes(coef_main, .data[[paste0("coef_", analysis)]],
                   color = Exposure, shape = label)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.6) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
    comparison_style() +
    scale_shape_manual(
      values = c("Suggestive" = 15, "Bonferroni" = 17),
      breaks = c("Suggestive", "Bonferroni"), drop = FALSE
    ) +
    labs(x = "Beta in Main Model", y = "Beta in Sensitivity Model",
         color = "Exposure", shape = "CpG significance") +
    theme(panel.border = element_blank())
}

# Main vs sensitivity models ----------------------------------------------
Main_model <- read_effects("main")
save(Main_model, file = file.path(output_dir, "Main_model_compare_plot.RData"))

for (analysis in c("90th", "99th", "ap", "complications", "EUR", "celltype", "gwase")) {
  Sens <- read_effects(analysis)
  Main_Sens <- left_join(Main_model, Sens, by = c("CpG", "Exposure"))
  plot <- plot_sensitivity(Main_Sens, analysis)
  ggsave(file.path(plot_dir, paste0("main_", analysis, "_plot.png")),
         plot = plot, width = 12, height = 8, dpi = 600)
}

# Males vs females ---------------------------------------------------------
ewas_sig <- read.csv(file.path(annotation_dir, "Sig_DMPs_annotated.csv"))
ewas_sig <- ewas_sig %>%
  subset(Bonferroni == "yes" | Suggestive == "yes") %>%
  mutate(label = ifelse(Bonferroni == "yes", "Bonferroni", "Suggestive"))

ewas_sig_sex <- ewas_sig %>%
  filter(grepl("[FM]$", Exposure)) %>%
  mutate(Sex = ifelse(grepl("F$", Exposure), "Female", "Male")) %>%
  select(CpG, Sex, Exposure, label)
ewas_sig_sex$Exposure <- sub(" [FM]$", "", ewas_sig_sex$Exposure)

read_sex_effects <- function(sex) {
  rows <- lapply(seq_len(nrow(exposures)), function(i) {
    data <- read_qc(i, sex = sex) %>%
      select(CpG, coef, se, pval) %>%
      mutate(Exposure = exposures$Exposure[i]) %>%
      right_join(ewas_sig_sex, by = c("CpG", "Exposure"))
    rename_effects(data, sex)
  })
  do.call(rbind, rows) %>% filter(!is.na(.data[[paste0("coef_", sex)]]))
}

male <- read_sex_effects("M")
female <- read_sex_effects("F")
sex_BN <- male %>%
  left_join(female, by = c("CpG", "Exposure", "Sex", "label")) %>%
  mutate(coef_F = coef_F * 10, coef_M = coef_M * 10)

sex <- ggplot(sex_BN, aes(coef_M, coef_F, color = Exposure, shape = label)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.6) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  comparison_style() +
  scale_shape_manual(values = c("Suggestive" = 15, "Bonferroni" = 17), drop = FALSE) +
  coord_equal() +
  labs(x = expression(beta ~ "in males"), y = expression(beta ~ "in females"),
       color = "Exposure", shape = "CpG significance")
ggsave(file.path(plot_dir, "sex_plot.png"), plot = sex, width = 12, height = 8, dpi = 600)
