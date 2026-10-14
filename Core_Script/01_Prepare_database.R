#-----------------------------------------------------------------------------#
#                    prepare the exposure and covariates                      #
#-----------------------------------------------------------------------------#
######################
# 1) Load packages
######################
library(tidyverse)
library(readr)
library(dplyr)
library(splines)

setwd("~/010_EWAS_temp_YZ") # change to your path
######################
# 2) combine database
######################
###############         covariates
load("/db/Data_Prepare/covariates_neuro_imputed.RData") # prepare covariates for your study
###############         exposure
load("/db/heat_count_whole_pregnancy.RData") # prepare exposure for your study

Exposure_Covariates <- left_join(covariates_neuro_imputed, heat_count, by= "id_mother")
ns_mat <- ns(as.numeric(Exposure_Covariates$conception_date), df = 12)
colnames(ns_mat) <- paste0("ns_cd_", 1:ncol(ns_mat))
Exposure_Covariates <- cbind(Exposure_Covariates, ns_mat)

save(Exposure_Covariates, file = "/db/Exposure_Covariates.rds")

# combine with Basename database 
load("/db/db_DNAm_vars_25022025.RData")
phenodataframe <- left_join(db, Exposure_Covariates, by = c("id_mother","id_child")) # phenodata with the exposures and covariates

save(phenodataframe, file = "/db/phenodataframe.rds")
#################################################################################
# END