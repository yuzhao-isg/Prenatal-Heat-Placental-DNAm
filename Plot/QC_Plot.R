#-----------------------------------------------------------------------------#
#                                    QC Plots                                 #
#-----------------------------------------------------------------------------#
rm(list = ls())
######################
# 1) Load packages
######################
library(tidyverse)
library(EASIER)
library(dplyr)
setwd("~/010_EWAS_temp_YZ") # change to your path

######################
# 2) Load database
######################
result_files <- c(
  indoor_HI95 = "/results/Main_Model/CpG_quality_control/Indoor/Exclude_CpGs_results_indoor_HI95.RData",
  indoor_HI95_F = "/results/Main_Model/CpG_quality_control/Indoor/count_HI_exposure_95_indoor_day/Exclude_CpGs_results_indoor_HI95_F.RData",
  indoor_HI95_M = "/results/Main_Model/CpG_quality_control/Indoor/count_HI_exposure_95_indoor_day/Exclude_CpGs_results_indoor_HI95_M.RData",
  indoor_T95 = "/results/Main_Model/CpG_quality_control/Indoor/Exclude_CpGs_results_indoor_T95.RData",
  indoor_T95_F = "/results/Main_Model/CpG_quality_control/Indoor/count_heat_exposure_95_indoor_day/Exclude_CpGs_results_indoor_T95_F.RData",
  indoor_T95_M = "/results/Main_Model/CpG_quality_control/Indoor/count_heat_exposure_95_indoor_day/Exclude_CpGs_results_indoor_T95_M.RData",
  outdoor_HI95 = "/results/Main_Model/CpG_quality_control/Outdoor/Exclude_CpGs_results_outdoor_HI95.RData",
  outdoor_HI95_F = "/results/Main_Model/CpG_quality_control/Outdoor/count_HI_exposure_95_outdoor_day/Exclude_CpGs_results_outdoor_HI95_F.RData",
  outdoor_HI95_M = "/results/Main_Model/CpG_quality_control/Outdoor/count_HI_exposure_95_outdoor_day/Exclude_CpGs_results_outdoor_HI95_M.RData",
  outdoor_T95 = "/results/Main_Model/CpG_quality_control/Outdoor/Exclude_CpGs_results_outdoor_T95.RData",
  outdoor_T95_F = "/results/Main_Model/CpG_quality_control/Outdoor/count_heat_exposure_95_outdoor_day/Exclude_CpGs_results_outdoor_T95_F.RData",
  outdoor_T95_M = "/results/Main_Model/CpG_quality_control/Outdoor/count_heat_exposure_95_outdoor_day/Exclude_CpGs_results_outdoor_T95_M.RData"
)

for (analysis_name in names(result_files)) {
  load(result_files[[analysis_name]])
  assign(analysis_name, get(paste0("result_", analysis_name)) %>%
    mutate(Exposure = analysis_name))
}
rm(result_files, analysis_name)
################################################################################


######################
# 3) Distribution Plot
######################
plot_specs <- expand.grid(
  location = c("indoor", "outdoor"),
  exposure = c("HI95", "T95"),
  sex = c("", "F", "M"),
  stringsAsFactors = FALSE
)
plot_specs$object_name <- paste0(
  plot_specs$location, "_", plot_specs$exposure,
  ifelse(plot_specs$sex == "", "", paste0("_", plot_specs$sex))
)

save_qc_plot <- function(results, filename, plot_type, title) {
  png(filename, width = 2000, height = 1600, res = 300)
  on.exit(dev.off())

  if (plot_type == "distribution") {
    plot_distribution(results$pval, main = title, xlab = "P_value")
  } else {
    qqman::qq(
      results$pval,
      main = sprintf(
        "QQ-plot of %s (lambda = %f)",
        title,
        get_lambda(results, "pval")
      )
    )
  }
}

for (index in seq_len(nrow(plot_specs))) {
  plot_spec <- plot_specs[index, ]
  plot_results <- get(plot_spec$object_name)
  location_label <- ifelse(plot_spec$location == "indoor", "Indoor", "Outdoor")
  exposure_label <- ifelse(plot_spec$exposure == "HI95", "Heat Index", "Temperature")
  sex_label <- ifelse(plot_spec$sex == "F", "Female", ifelse(plot_spec$sex == "M", "Male", ""))
  plot_title <- paste(location_label, exposure_label)
  if (nzchar(sex_label)) {
    plot_title <- paste(plot_title, sex_label)
  }
  plot_suffix <- ifelse(
    plot_spec$sex == "F", "_Female",
    ifelse(plot_spec$sex == "M", "_Male", "")
  )
  output_directory <- file.path(
    "/results/Main_Model/QC_Plot",
    ifelse(plot_spec$location == "indoor", "Indoor", "Outdoor")
  )

  save_qc_plot(
    plot_results,
    file.path(output_directory, paste0("Pvalue_", plot_spec$location, "_", plot_spec$exposure, plot_suffix, ".png")),
    "distribution",
    paste("Pvalue of", plot_title)
  )
  save_qc_plot(
    plot_results,
    file.path(output_directory, paste0("QQ_", plot_spec$location, "_", plot_spec$exposure, plot_suffix, ".png")),
    "qq",
    plot_title
  )
}

rm(plot_specs, save_qc_plot, index, plot_spec, plot_results, location_label,
   exposure_label, sex_label, plot_title, plot_suffix, output_directory)