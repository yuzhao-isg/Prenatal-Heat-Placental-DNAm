#-----------------------------------------------------------------------------#
#                           Enrichment - Heatmap                              #
#-----------------------------------------------------------------------------#
# 1. CromStates_Placenta                                                      #
# 2. PMD_Placenta                                                             #
#-----------------------------------------------------------------------------# 

######################
# Load packages
######################
library(dplyr)
library(tidyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(ggplot2)


######################
# 1) CromStates_Placenta
######################
root_dir <- "/results/Downstream/Enrichment/CromStates"
#-------------------------------------------------------
# 1. Sorted according to the Roadmap 15-state model
#    1–8: Active
#    9–15: Inactive/repressed
#-------------------------------------------------------

CromStates_enrich_results <- read.csv(file.path(root_dir, "CromStates_enrich_results.csv"), check.names = FALSE)
state_map <- c(
  "1_TssA" = "Active TSS",
  "2_TssAFlnk" = "Flanking active TSS",
  "10_TssBiv" = "Bivalent TSS",
  "7_Enh" = "Active enhancer",
  "6_EnhG" = "Genic enhancer",
  "12_EnhBiv" = "Bivalent enhancer",
  "4_Tx" = "Transcribed region",
  "3_TxFlnk" = "Transcription flanking region",
  "5_TxWk" = "Weak transcription",
  "13_ReprPC" = "Repressed polycomb",
  "14_ReprPCWk" = "Weak repressed Polycomb",
  "8_ZNF/Rpts" = "ZNF genes & repeats",
  "9_Het" = "Heterochromatin",
  "15_Quies" = "Quiescent",
  "11_BivFlnk" = "Bivalent Flanking"
)
state_order <- c(
  # Active states
  "Active TSS",
  "Flanking active TSS",
  "Transcription flanking region",
  "Transcribed region",
  "Weak transcription",
  "Genic enhancer",
  "Active enhancer",
  "ZNF genes & repeats",
  
  # Inactive/repressed states
  "Heterochromatin",
  "Bivalent TSS",
  "Bivalent Flanking",
  "Bivalent enhancer",
  "Repressed polycomb",
  "Weak repressed Polycomb",
  "Quiescent"
)

state_activity_map <- setNames(
  c(
    rep("Active", 8),
    rep("Inactive/repressed", 7)
  ),
  state_order
)
activity_levels <- c(
  "Active",
  "Inactive/repressed"
)
activity_colors <- c(
  "Active" = "#009E73",
  "Inactive/repressed" = "#7A7A7A"
)

#------------------------------
# 2. Add chromatin-state names
#------------------------------

df <- CromStates_enrich_results %>%
  mutate(
    chrom_state = recode(Data, !!!state_map),
    chrom_state = factor(chrom_state, levels = state_order)
  )

#------------------------------
# 3. Unified plotting function
#-------------------------------

make_chromstate_heatmap <- function(data,
                                    analysis_type,
                                    plot_title,
                                    output_file) {
  
  plot_df <- data %>%
    filter(Expression_level == analysis_type) %>%
    mutate(
      OR_sig = if_else(
        !is.na(FDR_p) & FDR_p < 0.05,
        log2(OR),
        NA_real_
      ),
      panel = paste(exposure, group, sep = "_"),
      chrom_state = as.character(chrom_state)
    ) %>%
    select(panel, chrom_state, OR_sig)
  
  # Check for duplicate results in the same heatmap cell
  duplicated_cells <- plot_df %>%
    count(panel, chrom_state) %>%
    filter(n > 1)
  
  if (nrow(duplicated_cells) > 0) {
    stop(
      "Some panel × chromatin-state combinations are duplicated. ",
      "Please inspect duplicated_cells."
    )
  }
  
  panel_order <- unique(plot_df$panel)
  plot_df_complete <- expand_grid(
    chrom_state = state_order,
    panel = panel_order
  ) %>%
    left_join(
      plot_df,
      by = c("chrom_state", "panel")
    ) %>%
    mutate(
      chrom_state = factor(
        chrom_state,
        levels = state_order
      ),
      panel = factor(
        panel,
        levels = panel_order
      )
    ) %>%
    arrange(chrom_state, panel)
  
  # Convert to matrix
  mat <- plot_df_complete %>%
    pivot_wider(
      names_from = panel,
      values_from = OR_sig
    ) %>%
    column_to_rownames("chrom_state") %>%
    as.matrix()
  
  # The active/inactive classification corresponding to each row
  state_activity <- factor(
    unname(state_activity_map[rownames(mat)]),
    levels = activity_levels
  )
  
  # Left-side color legend
  activity_annotation <- rowAnnotation(
    `Chromatin-state group` = state_activity,
    col = list(
      `Chromatin-state group` = activity_colors
    ),
    simple_anno_size = unit(5, "mm"),
    show_annotation_name = FALSE,
    annotation_legend_param = list(
      `Chromatin-state group` = list(
        title = "Chromatin-state group"
      )
    )
  )
  breaks <- seq(-3, 3, by = 0.25)
  
  heatmap_colors <- colorRampPalette(
    c("#2166AC", "white", "#B2182B")
  )(length(breaks))
  
  col_fun <- colorRamp2(
    breaks,
    heatmap_colors
  )
  
  ht <- Heatmap(
    mat,
    name = "log2(Odds Ratio)",
    col = col_fun,
    na_col = "white",
    
    # Active/inactive 
    left_annotation = activity_annotation,
    row_split = state_activity,
    row_gap = unit(2.5, "mm"),
    row_title = NULL,
    
    cluster_rows = FALSE,
    cluster_row_slices = FALSE,
    cluster_columns = FALSE,
    
    column_names_rot = 45,
    rect_gp = gpar(
      col = "grey80",
      lwd = 0.5
    ),
    
    column_title = plot_title,
    column_title_gp = gpar(
      fontsize = 14,
      fontface = "bold"
    ),
    
    width = unit(12, "cm"),
    height = unit(8, "cm")
  )
  
  pdf(
    output_file,
    width = 14,
    height = 10
  )
  
  draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    padding = unit(c(8, 8, 25, 25), "mm")
  )
  
  dev.off()
  
  invisible(ht)
}
##################

hyper <- make_chromstate_heatmap(
  data = df,
  analysis_type = "Fisher_Hyper",
  plot_title = paste0(
    "Chromatin-state enrichment of hypermethylated CpGs ",
    "(FDR < 0.05)"
  ),
  output_file = file.path(
    root_dir,
    "heatmap_hyper_grouped.pdf"
  )
)

hypo <- make_chromstate_heatmap(
  data = df,
  analysis_type = "Fisher_Hypo",
  plot_title = paste0(
    "Chromatin-state enrichment of hypomethylated CpGs ",
    "(FDR < 0.05)"
  ),
  output_file = file.path(
    root_dir,
    "heatmap_hypo_grouped.pdf"
  )
)
#################################################################################


######################
# 2) PMD_Placenta
######################
root_dir <- "/results/Downstream/Enrichment/PMD_Pla"

PMD_enrich_results <- read.csv(file.path(root_dir, "PMD_enrich_results.csv"), check.names = FALSE)

PMD_enrich_results$group <- factor(
  PMD_enrich_results$group,
  levels = c("Overall", "Female", "Male")
)

PMD_enrich_results$Expression_level <-
  sub("^Fisher_", "", PMD_enrich_results$Expression_level)

PMD_enrich_results$exposure <- factor(
  PMD_enrich_results$exposure,
  levels = c(
    "Indoor Temperature",
    "Indoor Heat Index",
    "Outdoor Temperature",
    "Outdoor Heat Index"
  )
)

PMD_enrich_results <- PMD_enrich_results %>%
  filter(Data == "NotNA")

p <- ggplot(
  PMD_enrich_results,
  aes(
    y = exposure,
    x = OR,
    colour = Expression_level
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = 2,
    colour = "grey60"
  ) +
  geom_errorbarh(
    aes(xmin = OR.inf, xmax = OR.sup),
    position = position_dodge(width = 0.6),
    height = 0.15,
    linewidth = 0.8
  ) +
  geom_point(
    position = position_dodge(width = 0.6),
    size = 3
  ) +
  facet_grid(group ~ .) +
  scale_colour_manual(values = c(
    Hyper = "#D55E00",
    Hypo = "#0072B2"
  )) +
  scale_x_log10() +
  theme_bw(base_size = 12) +
  labs(
    x = "Odds ratio (95% CI)",
    y = NULL,
    colour = NULL
  ) +
  geom_text(
    aes(label = ifelse(FDR_p < 0.05, "*", "")),
    position = position_dodge(width = 0.6),
    vjust = -0.5,
    show.legend = FALSE
  ) +
  theme(
    panel.grid = element_blank(),
    legend.position = "top"
  )

ggsave(
  file.path(root_dir, "PMD_enrich_results.png"),
  plot = p,
  width = 8,
  height = 6,
  dpi = 600,
  limitsize = FALSE
)
################################################################################
# END