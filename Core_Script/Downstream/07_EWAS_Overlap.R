#-----------------------------------------------------------------------------#
#                          Downstream - EWAS Overlap                          #
#-----------------------------------------------------------------------------#
# 1. EWAS Catalog                                                             #
# 2. EWAS Atlas                                                               #
#-----------------------------------------------------------------------------# 

rm(list = ls())
setwd("~/010_EWAS_temp_YZ") # change to your path
######################
# Load packages
######################
# install.packages("devtools")
#devtools::install_github("MRCIEU/ewascatalog-r")
library(readr)
library(dplyr)
library(readxl)
library(data.table)
library(ewascatalog)

######################
# 1) EWAS catalog - Run in HPC
######################
ewas_sig <- read.csv("/results/Downstream/Annotation/Sig_DMPs_annotated.csv")
sig_cpgs<-ewas_sig$CpG
sig_cpgs<-sig_cpgs[!duplicated(sig_cpgs)]

# Download the ewas catalog-results
download.file( url = "https://www.ewascatalog.org/static//docs/ewascatalog-results.txt.gz", 
               destfile = "/results/Downstream/EWAS_Overlap/EWAS_catalogue/ewascatalog_results.txt.gz", 
               method = "curl", extra = "-k")
# Download the ewas catalog-meta-data
# https://www.ewascatalog.org/download/

##### load results data
ewascat_db <- fread("/results/Downstream/EWAS_Overlap/EWAS_catalogue/ewascatalog_results.txt.gz")
colnames(ewascat_db)
cat("number of row:", nrow(ewascat_db), "\n")
cat("column name:", colnames(ewascat_db), "\n")
# there are two column name cpg
colnames(ewascat_db) <- c( "cpg", "beta", "se", "p", "details", "study_id", "cpg_dup", "chrpos", "chr", "pos", "gene", "type", "assocs" )
ewascat_results <- ewascat_db[, .(cpg, beta, se, p, study_id, chr, pos, gene)]
overlap <- ewascat_db[ewascat_results$cpg %in% sig_cpgs, ]
write.csv(overlap, "/results/Downstream/EWAS_Overlap/EWAS_catalogue/EWAS_catalog_overlap.csv", row.names = FALSE)

##### load meta-data
ewascat_studies <- fread("/results/Downstream/EWAS_Overlap/EWAS_catalogue/ewascatalog-studies.txt.gz")
# combine two database by study_id
ewascat_final <- overlap %>%
  left_join(ewascat_studies, by = "study_id")
write.csv(ewascat_final, "/results/Downstream/EWAS_Overlap/EWAS_catalogue/EWAS_catalog_final.csv", row.names = FALSE)
################################################################################


######################
# 2) EWAS Atlas 
######################

ewas_sig <- read.csv("/results/Downstream/Annotation/Sig_DMPs_annotated.csv")
sig_cpgs<-ewas_sig$CpG
sig_cpgs<-sig_cpgs[!duplicated(sig_cpgs)]
write.table(sig_cpgs,
            file = "/results/Downstream/EWAS_Overlap/EWAS_Atlas/sig_cpgs.txt",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

# submit the txt file to:https://ngdc.cncb.ac.cn/ewas/atlas
# https://ngdc.cncb.ac.cn/ewas/toolkit
################################################################################
# END