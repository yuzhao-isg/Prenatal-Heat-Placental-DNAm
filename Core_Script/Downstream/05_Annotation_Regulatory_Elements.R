#-----------------------------------------------------------------------------#
#               Downstream - Annotation Regulatory Elements                   #
#-----------------------------------------------------------------------------#
# 1. placental cis methylation quantitative trait loci (cis-mQTLs)            #
# 2. partially methylated domains (PMDs)                                      #
# 3. germline differentially methylated regions (gDMRs)                       #
# 4. Chromatin states                                                         #
# 5. imprinting control regions (ICRs)                                        #
#-----------------------------------------------------------------------------# 
rm(list = ls())
setwd("~/010_EWAS_temp_YZ") # change to your path
# This script should better run in HPC
######################
# Load packages
######################
library(dplyr)
library(data.table)
library(readr)
library(GenomicRanges)
library(tidyr)
library(stringr)

######################
# Load database
######################
ewas_sig <- read.csv("/results/Downstream/Annotation/Sig_DMPs_annotated.csv")
ewas_sig <- ewas_sig %>%
  select(CpG,Exposure,coef,se,pval,bonferroni_pval, Bonferroni, FDR_pval, FDR, Suggestive) 
any(duplicated(ewas_sig_all[, c("CpG")]))
sum(duplicated(ewas_sig_all[, c("CpG")]))

######################
# 1) cis-mQTLs
######################
#- mQTLs de PACE 
# https://pace-placenta-mqtl.streamlit.app/Resources
mQTLs_PACE <- read_tsv("/db/cis_mQTLs/gwama_filt_all_nstudies_gt1_topsnp.tsv")
colnames(mQTLs_PACE)
# [1] "rs_number"        "cpg"              "snp"              "reference_allele" "other_allele"     "eaf"              "beta"             "se"              
# [9] "beta_95L"         "beta_95U"         "z"                "p-value"          "_-log10_p-value"  "q_statistic"      "q_p-value"        "i2"              
# [17] "n_studies"        "n_samples"        "effects"  
mQTLs_pace<-left_join(ewas_sig_all,mQTLs_PACE,by = c("CpG" = "cpg"))

mQTLs_by_CpG_pace <- mQTLs_pace %>%
  group_by(CpG, Exposure) %>%
  summarise(mQTL_pace = paste(snp, collapse = ";"))

write.table(mQTLs_by_CpG_pace,"/results/Downstream/Genomic_regulatory/mQTLs_by_CpG_pace.txt",sep="\t")


#------------------------------------------
# PMD, gDMR and Chromatin states
#------------------------------------------
################################################################
#------- EASIER FUNCTIONS-DO NOT EDIT--------------------------
getEnrichGenomicRanges <- function(seqnames, start, end)
{
  # Creates the name to assign to GRanges
  pmdname <- paste(seqnames, start, end, sep = ';')
  
  Enrich.GRange <- GRanges(
    seqnames = Rle(seqnames),
    ranges = IRanges(start, end=end),
    name = pmdname
  )
  names(Enrich.GRange) <- Enrich.GRange$name
  
  return(Enrich.GRange)
  
}

findOverlapValues <- function(selected, complete,  outputdir = ".", outputfile = NULL)
{
  # Find overlaps
  ranges <- subsetByOverlaps(selected, complete)
  hits <- findOverlaps(selected, complete)
  
  # Find the names of states associated to each position
  idx <- subjectHits(hits)
  values <- DataFrame( States  = complete$name[idx])
  qh <- queryHits(hits)
  
  ans <- list("ranges" = ranges,
              "hits" = hits,
              "shits" = idx,
              "qhits" = qh,
              "values" = values)
  
  return(ans)
  
}

#############################################################

######################
# 2) PMD
######################

snp$CpG_chr_numeric <- as.numeric(gsub("chr", "", snp$CpG_chrm))


# Convert to Genomic Ranges
data.GRange <- GRanges(
  seqnames = Rle(snp$CpG_chrm),
  ranges=IRanges(snp$CpG_beg, end=snp$CpG_end),
  name=snp$CpG,
  chr=snp$CpG_chrm,
  pos=snp$CpG_beg
)
names(data.GRange) <- data.GRange$name
# Create genomic ranges from PMD data
#Load db with PMD (https://github.com/isglobal-brge/EASIER/tree/main/data)
load("/db/PMD_placenta.RData")

PMD.GRange <- getEnrichGenomicRanges(PMD_placenta$Chr_PMD, PMD_placenta$Start_PMD, PMD_placenta$End_PMD)

# Find overlaps between CpGs and PMD (find subject hits, query hits )
overPMD <- findOverlapValues(data.GRange, PMD.GRange)

#Create a data.frame with CpGs and PMDs information
mdata <- cbind.data.frame(DataFrame(CpG = data.GRange$name[overPMD$qhits]), DataFrame(PMD = PMD.GRange$name[overPMD$shits]))
any(duplicated(mdata[, c("CpG", "PMD")]))
sum(duplicated(mdata[, c("CpG", "PMD")]))
mdata_unique <- mdata[!duplicated(mdata[, c("CpG", "PMD")]), ] # REMOVE duplicated rows
# Merge with results of mQTLs
mQTLs <- left_join(ewas_sig, mQTLs_by_CpG_pace, by = c("CpG","Exposure"))
crom_data <- left_join(mQTLs, mdata_unique, by = "CpG")
nrow(crom_data)
write.table(crom_data,"/results/Downstream/Genomic_regulatory/comp_mQTLs_pmd.txt",sep=",")


######################
# 3) gDMR
######################
comp_pmd<-crom_data

data.GRange <- GRanges(
  seqnames = Rle(comp_pmd$CpG_chrm),
  ranges=IRanges(comp_pmd$CpG_beg, end=comp_pmd$CpG_end),
  name=comp_pmd$CpG,
  chr=comp_pmd$CpG_chrm,
  pos=comp_pmd$CpG_beg
)
names(data.GRange) <- data.GRange$name

# Create genomic ranges from gDMR (imprinted regions) 

#load gDMR (https://github.com/isglobal-brge/EASIER/tree/main/data)

load("/db/IR_Placenta.RData")

# Create genomic ranges from DMR data
DMR.GRange <- getEnrichGenomicRanges(IR_Placenta$Chr_DMR, IR_Placenta$Start_DMR, IR_Placenta$End_DMR)

# Find overlaps between CpGs and DMR (find subject hits, query hits )
overDMR <- findOverlapValues(data.GRange, DMR.GRange )

#Create a data.frame with CpGs and DMRs information
mdata <- cbind.data.frame(DataFrame(CpG = data.GRange$name[overDMR$qhits]), DataFrame(gDMR = DMR.GRange$name[overDMR$shits]))
head(mdata)
any(duplicated(mdata[, c("CpG", "gDMR")]))
sum(duplicated(mdata[, c("CpG", "gDMR")]))
mdata_unique <- mdata[!duplicated(mdata[, c("CpG", "gDMR")]), ]

# Merge with results from meta-analysis (A2)
crom_data <- left_join(comp_pmd, mdata_unique, by = "CpG")
nrow(crom_data)

write.table(crom_data,"/results/Downstream/Genomic_regulatory/comp_mQTLs_pmd_gdmr.txt",sep=",")


######################
# 4) Chromatine states
######################

comp_pmd_gdmr<-crom_data

data.GRange <- GRanges(
  seqnames = Rle(comp_pmd_gdmr$CpG_chrm),
  ranges=IRanges(comp_pmd_gdmr$CpG_beg, end=comp_pmd_gdmr$CpG_end),
  name=comp_pmd_gdmr$CpG,
  chr=comp_pmd_gdmr$CpG_chrm,
  pos=comp_pmd_gdmr$CpG_beg
)
names(data.GRange) <- data.GRange$name

#load chromatin states (https://github.com/isglobal-brge/EASIER/tree/main/data)

load("/db/FP_15_E091.RData")
# Find overlaps between CpGs and Fetal Placenta (States 15 and 18)
over15 <- findOverlapValues(data.GRange, FP_15_E091 )

#Create a data.frame with CpGs and DMRs information
mdata <- cbind.data.frame(DataFrame(CpG = data.GRange$name[over15$qhits]), DataFrame(ChrmSt15 = FP_15_E091$name[over15$shits]))

head(mdata)
any(duplicated(mdata[, c("CpG", "ChrmSt15")]))
sum(duplicated(mdata[, c("CpG", "ChrmSt15")]))
mdata_unique <- mdata[!duplicated(mdata[, c("CpG", "ChrmSt15")]), ]

mdata_unique_by_CpG <- mdata_unique %>%
  group_by(CpG) %>%
  summarise(ChrmSt15 = paste(ChrmSt15, collapse = ";"))

# Merge with results from meta-analysis (A2)
crom_data <- left_join(comp_pmd_gdmr, mdata_unique_by_CpG, by = "CpG")
nrow(crom_data)

write.table(crom_data,"/results/Downstream/Genomic_regulatory/comp_mQTLs_pmd_gdmr_chrmst.txt",sep=",")


######################
# 5) imprinting regions
######################
# Input:
#   1. Jima et al., 2022: ICR-level annotation
#   2. CpG sites mapped to ICRs, including tissue annotations
#   3. GeneImprint: imprinted-gene annotation
#
# Main analysis: Jima 2022
# Carreras-Gallo 2024: retained for reference, not run

# ============================================================
# 5.1. Main analysis: Jima 2022 ICR annotation
# ============================================================

# ---- Read CpG-to-ICR mapping ----

jima_cpg <- read_excel(
  "db/imprinting_regions/Jima2022_CpG_ICR.xlsx"
) %>%
  mutate(
    CpG_chr = if_else(
      str_detect(as.character(CpG_chr), "^chr"),
      as.character(CpG_chr),
      paste0("chr", CpG_chr)
    ),
    CpG_start = as.integer(CpG_start)
  ) %>%
  rename(CpG_end_ICR = CpG_end)

# ---- Read Jima 2022 ICR-level annotation ----

jima_icr <- read_excel(
  "db/imprinting_regions/Jima2022_1488_ICRs.xlsx"
) %>%
  mutate(
    # Remove symbols used to flag ICR entries
    ID = str_remove_all(as.character(ID), "[*^#]"),
    
    # Extract ICR genomic coordinates
    Jima_chr = str_extract(`Genomic Coordinates`, "^chr[^:]+"),
    Jima_start = as.integer(
      str_match(
        `Genomic Coordinates`,
        "^chr[^:]+:(\\d+)-(\\d+)$"
      )[, 2]
    ),
    Jima_end = as.integer(
      str_match(
        `Genomic Coordinates`,
        "^chr[^:]+:(\\d+)-(\\d+)$"
      )[, 3]
    )
  ) %>%
  select(
    ID,
    Jima_chr,
    Jima_start,
    Jima_end,
    `Parental Origin of Methylation`,
    `Nearest Transcript`,
    `Distance to Nearest Transcript`
  ) %>%
  distinct()


# ---- Match EWAS-significant CpGs to ICRs ----

ewas_icr_jima <- ewas_sig %>%
  mutate(
    CpG_chrm = if_else(
      str_detect(as.character(CpG_chrm), "^chr"),
      as.character(CpG_chrm),
      paste0("chr", CpG_chrm)
    ),
    CpG_beg = as.integer(CpG_beg)
  ) %>%
  left_join(
    jima_cpg,
    by = c(
      "CpG_chrm" = "CpG_chr",
      "CpG_beg" = "CpG_start"
    )
  ) %>%
  left_join(
    jima_icr,
    by = c("ICR_id" = "ID")
  )

# ---- Export results ----

output_dir <- "results/Downstream/Genomic_regulatory"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  ewas_icr_jima,
  file.path(output_dir, "ICR_Jima2022.csv"),
  row.names = FALSE,
  na = ""
)


# ============================================================
# 5.2. Additional annotation: GeneImprint
# ============================================================

# Annotate EWAS-significant CpG-associated genes using
# the GeneImprint database.
# This is a gene-level annotation, not an ICR overlap analysis.

geneimprint <- read_excel(
  "db/imprinting_regions/Geneimprint_web.xlsx"
) %>%
  mutate(
    across(
      where(is.character),
      ~ .x %>%
        str_replace_all("\u00C2", "") %>%
        str_replace_all("\u00A0", " ") %>%
        str_squish() %>%
        na_if("")
    )
  ) %>%
  distinct(Gene, .keep_all = TRUE)


# ---- Match annotated genes ----

ewas_geneimprint <- ewas_sig %>%
  mutate(orig_row_id = row_number()) %>%
  separate_rows(UCSC_RefGene_Name, sep = ";") %>%
  mutate(
    UCSC_RefGene_Name = str_trim(UCSC_RefGene_Name),
    UCSC_RefGene_Name = na_if(UCSC_RefGene_Name, "")
  ) %>%
  left_join(
    geneimprint,
    by = c("UCSC_RefGene_Name" = "Gene")
  ) %>%
  mutate(
    Location = na_if(str_trim(Location), ""),
    Status = na_if(str_trim(Status), ""),
    imprinting_anno = case_when(
      !is.na(Status) & !is.na(Location) ~
        paste0(Status, " (", Location, ")"),
      !is.na(Status) ~ Status,
      !is.na(Location) ~ Location,
      TRUE ~ NA_character_
    )
  )


# ---- Restore one row per original EWAS result ----

ewas_cols <- names(ewas_sig)

ewas_geneimprint_final <- ewas_geneimprint %>%
  group_by(orig_row_id) %>%
  summarise(
    across(
      all_of(setdiff(ewas_cols, "UCSC_RefGene_Name")),
      first
    ),
    UCSC_RefGene_Name = paste(
      unique(na.omit(UCSC_RefGene_Name)),
      collapse = ";"
    ),
    imprinting_anno = paste(
      unique(na.omit(imprinting_anno)),
      collapse = ";"
    ),
    .groups = "drop"
  ) %>%
  mutate(
    UCSC_RefGene_Name = na_if(UCSC_RefGene_Name, ""),
    imprinting_anno = na_if(imprinting_anno, "")
  ) %>%
  select(-orig_row_id)


# ---- Export results ----

write.csv(
  ewas_geneimprint_final,
  file.path(output_dir, "ICR_GeneImprint.csv"),
  row.names = FALSE,
  na = ""
)


# ============================================================
# 5.3. Carreras-Gallo et al., 2024
# ============================================================

# NOT RUN:
# The Carreras-Gallo file is retained as an alternative
# CpG-to-ICR mapping resource for future reference.
#
# Before using it, compare its probe coverage, ICR definitions,
# genome build, and overlap with the Jima-based mapping.
#
# carreras_icr <- read_excel(
#   "db/imprinting_regions/Carreras-Gallo2024_9757_ICRs_to_start.xlsx"
# )
#
# This resource is not used in the current main analysis.
################################################################################
# END
