###############################################
# Enrichment with MissMethyl using EnrichR databases
# GSEA
# https://bioconductor.org/packages//release/bioc/vignettes/methylGSA/inst/doc/methylGSA-vignette.html
rm(list = ls())
setwd("~/010_EWAS_temp_YZ") # change to your path

library(data.table)
library(methylGSA)
library(tidyverse)
library(org.Hs.eg.db)
library(dplyr)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

######################
# 1) Load datasets
######################
dirs <- list.dirs("/results/Main_Model/CpG_quality_control", full.names = TRUE)
filenames <- list.files(dirs, pattern = "*95.RData", full.names = TRUE)

dataAll <- c()

for (r in filenames) {
  obj_name <- load(r)
  main <- get(obj_name) %>%
    dplyr::select(CpG, coef, se, pval, Bonferroni, FDR_pval, CpG_chrm, CpG_beg, CpG_end) %>%
    dplyr::rename(bn = Bonferroni, fdr = FDR_pval)
  rownames(main) <- main$CpG

  file_base <- tools::file_path_sans_ext(basename(r))
  parts <- unlist(strsplit(file_base, "_"))
  name <- paste(parts[(length(parts) - 1):length(parts)], collapse = "_")

  assign(name, main)
  dataAll <- c(dataAll, name)
}
rm(main)

######################
# 2) Prepare gene annotation
######################
annotdat <- as.data.frame(getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))
annotdat <- annotdat[, c("Name", "pos", "chr", "UCSC_RefGene_Name")]
annotdat$ArrayEPIC <- "Yes"
rownames(annotdat) <- annotdat$Name
anno <- annotdat
rm(annotdat)

anno <- anno %>%
  as.data.frame() %>%
  rownames_to_column(var = "ProbeID")

prepareAnno.input <- anno %>%
  dplyr::select(ProbeID, UCSC_RefGene_Name) %>%
  mutate(Gene = str_split_i(string = UCSC_RefGene_Name, pattern = ";", i = 1)) %>%
  dplyr::rename(CpG = ProbeID) %>%
  dplyr::select(!UCSC_RefGene_Name)

FullAnnot <- prepareAnnot(prepareAnno.input)

######################
# 3) Load gene-set databases
######################
# Enrichr website: https://maayanlab.cloud/Enrichr/#libraries
read_enrichr_database <- function(file_path) {
  dat <- scan(file_path, what = "", sep = "\n")
  ldat <- strsplit(dat, "\t")

  names <- vapply(ldat, function(x) x[1], character(1))
  ldat <- lapply(ldat, function(x) x[-1])
  ldat <- lapply(ldat, function(x) x[x != ""])
  names(ldat) <- names

  lapply(ldat, function(x) {
    suppressMessages(
      AnnotationDbi::select(
        org.Hs.eg.db,
        keys = x,
        columns = c("ENTREZID", "SYMBOL"),
        keytype = "SYMBOL"
      )$ENTREZID
    )
  })
}

wiki <- read_enrichr_database("db/pathways/WikiPathways_2024_Human.txt")
reactome <- read_enrichr_database("db/pathways/Reactome_2022.txt")
biocarta <- read_enrichr_database("db/pathways/BioCarta_2016.txt")
panther <- read_enrichr_database("db/pathways/Panther_2016.txt")
Msigdb <- read_enrichr_database("db/pathways/MSigDB_Hallmark_2020.txt")
tf <- read_enrichr_database("db/pathways/ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.txt")
disgenet <- read_enrichr_database("db/pathways/DisGeNET.txt")
dsigdb <- read_enrichr_database("db/pathways/DSigDB.txt")
omimexp <- read_enrichr_database("db/pathways/OMIM_Expanded.txt")
omimdis <- read_enrichr_database("db/pathways/OMIM_Disease.txt")
gwascat <- read_enrichr_database("db/pathways/GWAS_Catalog_2025.txt")

#-----------------------------------------------------------------------------#
#                                   Enrichment                                #
#-----------------------------------------------------------------------------#
######################################################
#·····················································

# loop for all datasets
for (z in 1:length(dataAll)) {
  r <- get(dataAll[z])
  cpg.pval <- r$pval
  names(cpg.pval) <- rownames(r)
  
  
  ##############################
  #############################
  #PATHWAYS methylGSA
  #############################
  #############################
  #Databases: GO, KEGG, WikiPathways, Reactome, BioCarta, Panter, and MSigDB HALLMARK
  
  #----
  #GO
  #----
  gst <- methylRRA(cpg.pval = cpg.pval,
                   method = "GSEA",
                   minsize = 20,
                   maxsize = 500,
                   FullAnnot = FullAnnot,
                   GS.type = "GO")
  go.order<-gst[order(gst$pvalue),]
  
  #-----
  #KEGG 
  #-----
  
  gst<-methylRRA(cpg.pval = cpg.pval,
                 method = "GSEA",
                 minsize = 20,
                 maxsize = 500,
                 FullAnnot = FullAnnot,
                 GS.type = "KEGG")
  kegg.order<-gst[order(gst$pvalue),]
  
  #-----
  #WikiPathway_2023_Human
  #-----
  
  ### Specify arguments accordingly (using same arguments as with gometh function)
  ldat_entrez <- wiki
  gst <- methylRRA(cpg.pval = cpg.pval,
                   method = "GSEA",
                   minsize = 20,
                   maxsize = 500,
                   FullAnnot = FullAnnot,
                   GS.list = ldat_entrez,
                   GS.idtype = "ENTREZID")
  wiki.order<-gst[order(gst$pvalue),]
  
  #-----
  #Reactome_2022
  #-----
  
  ### Specify arguments accordingly (using same arguments as with gometh function)
  # ldat_entrez <- reactome
  gst <- methylRRA(cpg.pval = cpg.pval,
                   method = "GSEA",
                   minsize = 20,
                   maxsize = 500,
                   FullAnnot = FullAnnot,
                   GS.type = "Reactome")
  reactome.order<-gst[order(gst$pvalue),]
  
  #-----
  #BioCarta_2016
  #-----
  
  ### Specify arguments accordingly (using same arguments as with gometh function)
  ldat_entrez <- biocarta
  gst <- methylRRA(cpg.pval = cpg.pval,
                   method = "GSEA",
                   minsize = 20,
                   maxsize = 500,
                   FullAnnot = FullAnnot,
                   GS.list = ldat_entrez,
                   GS.idtype = "ENTREZID")
  biocarta.order<-gst[order(gst$pvalue),]
  
  #-----
  #Panther_2016.txt
  #-----
  
  ### Specify arguments accordingly (using same arguments as with gometh function)
  ldat_entrez <- panther
  gst <- methylRRA(cpg.pval = cpg.pval,
                   method = "GSEA",
                   minsize = 20,
                   maxsize = 500,
                   FullAnnot = FullAnnot,
                   GS.list = ldat_entrez,
                   GS.idtype = "ENTREZID")
  panther.order<-gst[order(gst$pvalue),]
  
  #-----
  #MSigDB
  #-----
  
  ### Specify arguments accordingly (using same arguments as with gometh function)
  ldat_entrez <- Msigdb
  gst <- methylRRA(cpg.pval = cpg.pval,
                   method = "GSEA",
                   minsize = 20,
                   maxsize = 500,
                   FullAnnot = FullAnnot,
                   GS.list = ldat_entrez,
                   GS.idtype = "ENTREZID")
  msigdb.order<-gst[order(gst$pvalue),]
  
  #########################
  #########################
  #TRANSCRIPTION FACTORS
  #########################
  #########################
  
  #-----
  #ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.txt
  #-----
  
  ### Specify arguments accordingly (using same arguments as with gometh function)
  ldat_entrez <- tf
  gst <- methylRRA(cpg.pval = cpg.pval,
                   method = "GSEA",
                   minsize = 20,
                   maxsize = 500,
                   FullAnnot = FullAnnot,
                   GS.list = ldat_entrez,
                   GS.idtype = "ENTREZID")
  tf.order<-gst[order(gst$pvalue),]
  
  ########################
  ########################
  #DISEASES
  ########################
  ########################
  
  #-----
  #DisGeNET.txt
  #-----
  
  ### Specify arguments accordingly (using same arguments as with gometh function)
  ldat_entrez <- disgenet
  gst <- methylRRA(cpg.pval = cpg.pval,
                   method = "GSEA",
                   minsize = 20,
                   maxsize = 500,
                   FullAnnot = FullAnnot,
                   GS.list = ldat_entrez,
                   GS.idtype = "ENTREZID")
  disgenet.order<-gst[order(gst$pvalue),]
  
  #-----
  #DSigDB.txt
  #-----
  
  ### Specify arguments accordingly (using same arguments as with gometh function)
  ldat_entrez <- dsigdb
  gst <- methylRRA(cpg.pval = cpg.pval,
                   method = "GSEA",
                   minsize = 20,
                   maxsize = 500,
                   FullAnnot = FullAnnot,
                   GS.list = ldat_entrez,
                   GS.idtype = "ENTREZID")
  dsigdb.order<-gst[order(gst$pvalue),]
  
  #-----
  #OMIM_Expanded.txt
  #-----
  
  ### Specify arguments accordingly (using same arguments as with gometh function)
  ldat_entrez <- omimexp
  gst <- methylRRA(cpg.pval = cpg.pval,
                   method = "GSEA",
                   minsize = 20,
                   maxsize = 500,
                   FullAnnot = FullAnnot,
                   GS.list = ldat_entrez,
                   GS.idtype = "ENTREZID")
  omimexp.order<-gst[order(gst$pvalue),]
  
  #-----
  #OMIM_Disease
  #-----
  
  ### Specify arguments accordingly (using same arguments as with gometh function)
  ldat_entrez <- omimdis
  gst <- methylRRA(cpg.pval = cpg.pval,
                   method = "GSEA",
                   minsize = 20,
                   maxsize = 500,
                   FullAnnot = FullAnnot,
                   GS.list = ldat_entrez,
                   GS.idtype = "ENTREZID")
  omimdis.order<-gst[order(gst$pvalue),]
  
  #-----
  #GWAS catalogue
  #-----
  
  ### Specify arguments accordingly (using same arguments as with gometh function)
  ldat_entrez <- gwascat
  gst <- methylRRA(cpg.pval = cpg.pval,
                   method = "GSEA",
                   minsize = 20,
                   maxsize = 500,
                   FullAnnot = FullAnnot,
                   GS.list = ldat_entrez,
                   GS.idtype = "ENTREZID")
  gwascat.order<-gst[order(gst$pvalue),]
  
  
  
  #######################
  #######################
  # save the data as rdata and xlsx
  #######################
  #######################
  
  save(go.order, kegg.order, wiki.order, reactome.order, biocarta.order, panther.order, msigdb.order, tf.order,
       disgenet.order, dsigdb.order, omimexp.order, omimdis.order, gwascat.order,
       file = paste("results/Downstream/Enrichment/GSEA/Overall/", dataAll[z], "_all_cpgs.Rdata", sep = ""))
  # 
  # for (s in c("go.order", "kegg.order", "wiki.order", "reactome.order", "biocarta.order", "panther.order", "msigdb.order", "tf.order",
  #             "disgenet.order", "dsigdb.order", "omimexp.order", "omimdis.order", "gwascat.order")) {
  #   tmp2 <- get(s)
  #   tmp2 <- cbind(rownames(tmp2), tmp2)
  #   colnames(tmp2)[1] <- 'Pathway'
  #   assign(paste(s, sep = ""), tmp2)
  # }
  # list_of_datasets <- list(go.order, kegg.order, wiki.order, reactome.order, biocarta.order, panther.order, msigdb.order, tf.order,
  #                          disgenet.order, dsigdb.order, omimexp.order, omimdis.order, gwascat.order)
  # names(list_of_datasets) <- c("GO", "KEGG", "WikiPathways", "Reactome", "BioCarta", "Panther", "MSigDB", "Transcription_Factors",
  #                              "DisGeNET", "DSigDB", "OMIM_Expanded", "OMIM_Diseases", "GWAS_Catalogue")
  # openxlsx::write.xlsx(list_of_datasets, file = paste("results/pathways/methylgsa", dataAll[z], ".xlsx", sep = ""))
}
