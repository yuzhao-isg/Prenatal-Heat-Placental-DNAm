###############################################
# Enrichment with MissMethyl using EnrichR databases
rm(list = ls())
setwd("~/010_EWAS_temp_YZ") # change to your path

library(data.table)
library(missMethyl)
library(org.Hs.eg.db)
library(dplyr)

######################
# 1) Loading results of interest
######################
dirs <- list.dirs("/results/Main_Model/CpG_quality_control", full.names = TRUE)
filenames <- list.files(dirs, pattern = "*95.RData", full.names = TRUE)

dataAll <- c()
for (r in filenames) {
  obj_name <- load(r)
  main <- get(obj_name) %>%
    select(CpG, coef, se, pval, Bonferroni, FDR_pval, CpG_chrm, CpG_beg, CpG_end) %>%
    rename(bn = Bonferroni, fdr = FDR_pval)
  rownames(main) <- main$CpG

  file_base <- tools::file_path_sans_ext(basename(r))
  parts <- unlist(strsplit(file_base, "_"))
  name <- paste(parts[(length(parts) - 1):length(parts)], collapse = "_")

  assign(name, main)
  dataAll <- c(dataAll, name)
}
rm(main)

######################
# 2) Load gene-set databases
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
library(msigdbr)

# loop for all datasets
for (name in dataAll[1:4]) {
  r <- get(name)
  
  #######################
  #######################
  # perform the analysis on suggestive CpGs
  #######################
  #######################
  
  tryCatch({cpgs <- getMappedEntrezIDs(r[which(r$pval < 1e-04), 'CpG'], r$CpG, array.type="EPIC")}, 
           error = function(e) { 
             cat("ERROR :",conditionMessage(e), "\n")
             cpgs <<- c()})
  if(length(cpgs) > 0) { 
    
    ##############################
    #############################
    #PATHWAYS MissMethyl
    #############################
    #############################
    #Databases: GO, KEGG, WikiPathways, Reactome, BioCarta, Panter, and MSigDB HALLMARK
    
    #----
    #GO
    #----
    
    gst<-gometh(r[which(r$pval < 1e-04), 'CpG'], r$CpG, collection = "GO", array.type="EPIC", sig.genes = T)
    go.order<-gst[order(gst$P.DE),]
    
    #-----
    #KEGG 
    #-----
    
    gst<-gometh(r[which(r$pval < 1e-04), 'CpG'], r$CpG, collection = "KEGG", array.type="EPIC", sig.genes = T)
    kegg.order<-gst[order(gst$P.DE),]
    
    #-----
    #WikiPathway_2023_Human
    #-----
    
    ### Specify arguments accordingly (using same arguments as with gometh function)
    ldat_entrez <- wiki
    gst <- gsameth(sig.cpg=r[which(r$pval < 1e-04), 'CpG'], all.cpg=r$CpG, collection=ldat_entrez, 
                   array.type="EPIC", plot.bias=FALSE, prior.prob=TRUE, anno=NULL, equiv.cpg=TRUE,
                   fract.counts=TRUE, genomic.features = "ALL", sig.genes = T)
    wiki.order<-gst[order(gst$P.DE),]
    
    #-----
    #Reactome_2022
    #-----
    
    ### Specify arguments accordingly (using same arguments as with gometh function)
    ldat_entrez <- reactome
    gst <- gsameth(sig.cpg=r[which(r$pval < 1e-04), 'CpG'], all.cpg=r$CpG, collection=ldat_entrez, 
                   array.type="EPIC", plot.bias=FALSE, prior.prob=TRUE, anno=NULL, equiv.cpg=TRUE,
                   fract.counts=TRUE, genomic.features = "ALL", sig.genes = T)
    reactome.order<-gst[order(gst$P.DE),]
    
    #-----
    #BioCarta_2016
    #-----
    
    ### Specify arguments accordingly (using same arguments as with gometh function)
    ldat_entrez <- biocarta
    gst <- gsameth(sig.cpg=r[which(r$pval < 1e-04), 'CpG'], all.cpg=r$CpG, collection=ldat_entrez, 
                   array.type="EPIC", plot.bias=FALSE, prior.prob=TRUE, anno=NULL, equiv.cpg=TRUE,
                   fract.counts=TRUE, genomic.features = "ALL", sig.genes = T)
    biocarta.order<-gst[order(gst$P.DE),]
    
    #-----
    #Panther_2016.txt
    #-----
    
    ### Specify arguments accordingly (using same arguments as with gometh function)
    ldat_entrez <- panther
    gst <- gsameth(sig.cpg=r[which(r$pval < 1e-04), 'CpG'], all.cpg=r$CpG, collection=ldat_entrez, 
                   array.type="EPIC", plot.bias=FALSE, prior.prob=TRUE, anno=NULL, equiv.cpg=TRUE,
                   fract.counts=TRUE, genomic.features = "ALL", sig.genes = T)
    panther.order<-gst[order(gst$P.DE),]
    
    #-----
    #MSigDB
    #-----
    
    ### Specify arguments accordingly (using same arguments as with gometh function)
    ldat_entrez <- Msigdb
    gst <- gsameth(sig.cpg=r[which(r$pval < 1e-04), 'CpG'], all.cpg=r$CpG, collection=ldat_entrez, 
                   array.type="EPIC", plot.bias=FALSE, prior.prob=TRUE, anno=NULL, equiv.cpg=TRUE,
                   fract.counts=TRUE, genomic.features = "ALL", sig.genes = T)
    msigdb.order<-gst[order(gst$P.DE),]
    
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
    gst <- gsameth(sig.cpg=r[which(r$pval < 1e-04), 'CpG'], all.cpg=r$CpG, collection=ldat_entrez, 
                   array.type="EPIC", plot.bias=FALSE, prior.prob=TRUE, anno=NULL, equiv.cpg=TRUE,
                   fract.counts=TRUE, genomic.features = "ALL", sig.genes = T)
    tf.order<-gst[order(gst$P.DE),]
    
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
    gst <- gsameth(sig.cpg=r[which(r$pval < 1e-04), 'CpG'], all.cpg=r$CpG, collection=ldat_entrez, 
                   array.type="EPIC", plot.bias=FALSE, prior.prob=TRUE, anno=NULL, equiv.cpg=TRUE,
                   fract.counts=TRUE, genomic.features = "ALL", sig.genes = T)
    disgenet.order<-gst[order(gst$P.DE),]
    
    #-----
    #DSigDB.txt
    #-----
    
    ### Specify arguments accordingly (using same arguments as with gometh function)
    ldat_entrez <- dsigdb
    gst <- gsameth(sig.cpg=r[which(r$pval < 1e-04), 'CpG'], all.cpg=r$CpG, collection=ldat_entrez, 
                   array.type="EPIC", plot.bias=FALSE, prior.prob=TRUE, anno=NULL, equiv.cpg=TRUE,
                   fract.counts=TRUE, genomic.features = "ALL", sig.genes = T)
    dsigdb.order<-gst[order(gst$P.DE),]
    
    #-----
    #OMIM_Expanded.txt
    #-----
    
    ### Specify arguments accordingly (using same arguments as with gometh function)
    ldat_entrez <- omimexp
    gst <- gsameth(sig.cpg=r[which(r$pval < 1e-04), 'CpG'], all.cpg=r$CpG, collection=ldat_entrez, 
                   array.type="EPIC", plot.bias=FALSE, prior.prob=TRUE, anno=NULL, equiv.cpg=TRUE,
                   fract.counts=TRUE, genomic.features = "ALL", sig.genes = T)
    omimexp.order<-gst[order(gst$P.DE),]
    
    #-----
    #OMIM_Disease
    #-----
    
    ### Specify arguments accordingly (using same arguments as with gometh function)
    ldat_entrez <- omimdis
    gst <- gsameth(sig.cpg=r[which(r$pval < 1e-04), 'CpG'], all.cpg=r$CpG, collection=ldat_entrez, 
                   array.type="EPIC", plot.bias=FALSE, prior.prob=TRUE, anno=NULL, equiv.cpg=TRUE,
                   fract.counts=TRUE, genomic.features = "ALL", sig.genes = T)
    omimdis.order<-gst[order(gst$P.DE),]
    
    #-----
    #GWAS catalogue
    #-----
    
    ### Specify arguments accordingly (using same arguments as with gometh function)
    ldat_entrez <- gwascat
    gst <- gsameth(sig.cpg=r[which(r$pval < 1e-04), 'CpG'], all.cpg=r$CpG, collection=ldat_entrez, 
                   array.type="EPIC", plot.bias=FALSE, prior.prob=TRUE, anno=NULL, equiv.cpg=TRUE,
                   fract.counts=TRUE, genomic.features = "ALL", sig.genes = T)
    gwascat.order<-gst[order(gst$P.DE),]
    
    
    
    #######################
    #######################
    # save the data as rdata 
    #######################
    #######################
    
    save(go.order, kegg.order, wiki.order, reactome.order, biocarta.order, panther.order, msigdb.order, tf.order, 
         disgenet.order, dsigdb.order, omimexp.order, omimdis.order, gwascat.order, 
         file = paste("results/Heat_days/Downstream/Enrichment/ORA/Overall/", name, "_e04.Rdata", sep = ""))
  }
}

# save xlsx
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
# openxlsx::write.xlsx(list_of_datasets, file = paste("results/Heat_days/Main_Results/Enrichment/pathways_all/", name, "_sugg_cpgs.xlsx", sep = ""))
