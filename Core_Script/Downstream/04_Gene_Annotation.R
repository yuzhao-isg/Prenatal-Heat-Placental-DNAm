#-----------------------------------------------------------------------------#
#                          Downstream - Gene Annotation                       #
#-----------------------------------------------------------------------------#
# 1. Positional genes  (Illumina)                                             #
# 2. Gene-enhancer in trophoblasts (EPIraction)                               #
# 3. Placental cis-eQTMs (PACE catalog)                                       #
# 4. Annotation for no Gene CpGs (not include in this study)                  #
#-----------------------------------------------------------------------------# 

rm(list = ls())
setwd("~/010_EWAS_temp_YZ") # change to your path
######################
# 1) Load packages
######################
library(dplyr)
library(EASIER)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(tibble)
library(readr)
library(data.table)
library(biomaRt)

######################
# 2) Load dataset
######################
result_files <- c(
  indoor_HI95 = "/results/Main_Model/CpG_quality_control/Indoor/Exclude_CpGs_results_indoor_HI95.RData",
  indoor_HI95_F = "/results/Main_Model/CpG_quality_control/Indoor/count_HI_exposure_95_indoor_day/Exclude_CpGs_results_indoor_HI95_F.RData",
  indoor_HI95_M = "/results/Main_Model/CpG_quality_control/Indoor/count_HI_exposure_95_indoor_day/Exclude_CpGs_results_indoor_HI95_M.RData",
  indoor_T95 = "/results/Main_Model/CpG_quality_control/Indoor/Exclude_CpGs_results_indoor_T95.RData",
  indoor_T95_F = "/results/Main_Model/CpG_quality_control/Indoor/count_heat_exposure_95_indoor_day/Exclude_CpGs_results_indoor_T95_F.RData",
  indoor_T95_M = "/results/Main_Model/CpG_quality_control/Indoor/count_heat_exposure_95_indoor_day/Exclude_CpGs_results_indoor_T95_M.RData",
  outdoor_HI95 = "/results/Main_Model/CpG_quality_control/Outdoor/Exclude_CpGs_results_outdoor_HI95.RData",
  outdoor_HI95_F = "/results/Main_Model/CpG_quality_control/Outdoor/count_HI_exposure_95_home/Exclude_CpGs_results_outdoor_HI95_F.RData",
  outdoor_HI95_M = "/results/Main_Model/CpG_quality_control/Outdoor/count_HI_exposure_95_home/Exclude_CpGs_results_outdoor_HI95_M.RData",
  outdoor_T95 = "/results/Main_Model/CpG_quality_control/Outdoor/Exclude_CpGs_results_outdoor_T95.RData",
  outdoor_T95_F = "/results/Main_Model/CpG_quality_control/Outdoor/count_heat_exposure_95_home/Exclude_CpGs_results_outdoor_T95_F.RData",
  outdoor_T95_M = "/results/Main_Model/CpG_quality_control/Outdoor/count_heat_exposure_95_home/Exclude_CpGs_results_outdoor_T95_M.RData"
)

for (analysis_name in names(result_files)) {
  load(result_files[[analysis_name]])
  assign(
    analysis_name,
    get(paste0("result_", analysis_name)) %>%
      mutate(Exposure = analysis_name)
  )
}
db_list <- mget(names(result_files))
rm(result_files, analysis_name)

######################
# 3) Combine the significant CpGs
######################
Sig_CpGs <- bind_rows(
  lapply(db_list, function(df) {
    df %>% filter(Bonferroni == "yes" | FDR == "yes" | Suggestive == "yes")
  })
)

Sig_CpGs <- Sig_CpGs %>%
  group_by(CpG) %>%
  mutate(CpG_repeat = ifelse(n() > 1, "repeat", "no")) %>%
  ungroup()

################################################################################

######################
# 4) Annotate a List of CpGs with Illumina 450K or EPIC data
######################
#' Get CpG annotations from Illumina 450K or EPIC array type
#'
#' @param CpGs string. Path with results from GWAMA
#' @param artype string, Illumina array type, 450K or EPIC
#' @param filename string vector with file name where annotations
#' @param outdir string vector. Route to output folder to store results
#'
#' @return dataframe with annotated CpGs
#'
#' @export
get_annotattions <- function(CpGs, artype, filename, outdir)
{
  
  # Get library for 450K or Epic
  if( toupper(artype) == '450K')
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  else
    library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  
  # Output filename
  outfilename <- tools::file_path_sans_ext(basename(filename))
  outfilename <- paste(outdir, outfilename, sep="/")
  
  # Download data
  if( toupper(artype) == '450K')
    ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  else
    ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  ann$gene <- rownames(ann)
  
  # Merge cpgs with annotations
  CpGs.annot <-  merge( as.data.frame(CpGs), ann, by.x ="CpG", by.y = "gene")
  
  # Remove data that we are not interested in (duplicate data) before write results
  ##..## CpGs.annot <- CpGs.annot[ ,-(1)]
  
  # Write data to a file
  write.table( CpGs.annot, paste0(outfilename,"_annotated.txt"), quote=F, sep="\t")
  
  # Return data
  return(CpGs.annot)
}

annotation_DMPs <- get_annotattions(Sig_CpGs, "EPIC", "Sig_DMPs", "/results/Main_Model/Downstream/Annotation/")
write.csv(annotation_DMPs, file = "/results/Main_Model/Downstream/Annotation/Sig_DMPs_annotated.csv", append = FALSE, col.names = TRUE)
###############################################################################


######################
# 5) Prepare txt file for EPIraction platform
######################
# https://epiraction.crg.eu/
Sig_CpGs <- read.delim("/results/Main_Model/Annotation/Sig_DMPs_annotated.txt", stringsAsFactors = FALSE)
cpg_positions <- Sig_CpGs[, c("chr", "pos", "CpG")]
con <- gzfile("/results/Main_Model/Downstream/Annotation/CpG_positions_hg38.txt.gz", "w")
writeLines("# CpG positions hg38", con)
write.table(cpg_positions,
            file = con,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

close(con)

# After using EPIraction platform to annotation
EPIraction<- read_tsv(
  "/results/Main_Model/Downstream/Annotation/EPIraction_annotated.txt",
  comment = "##", 
  col_types = cols()
)

######################
# 6) Enhancer from EPIraction platform 
######################
EPIraction <- read_tsv(
  "/results/Downstream/Annotation/EPIraction_annotated.txt",
  comment = "##", 
  col_types = cols()
)
# Only the EPIraction > 0.02 we consider as significant
EPIraction <- EPIraction %>%
  subset(EPIraction > 0.02)
any(duplicated(EPIraction[, c("#chrom", "position","annotation","enhancer","EPIraction","symbol")]))
sum(duplicated(EPIraction[, c("#chrom", "position","annotation","enhancer","EPIraction","symbol")]))
EPIraction_unique <- EPIraction[!duplicated(EPIraction[, c("#chrom", "position","annotation","enhancer","EPIraction","symbol")]), ]

EPIraction_final <- EPIraction_unique[!duplicated(EPIraction_unique[, c("#chrom", "position","annotation","enhancer","EPIraction","symbol")]), ]
save(EPIraction_final, file = "/results/Downstream/Annotation/EPIraction/EPIraction_final.RData")
################################################################################


######################
# 7) cis-eQTMs
######################
# load the data
ewas_sig <- read.csv("/results/Downstream/Annotation/Indoor/Sig_DMPs_annotated.csv")
ewas_sig_all <- ewas_sig %>%
  select(CpG,Exposure,coef,se,pval,bonferroni_pval, Bonferroni, FDR_pval, FDR, Suggestive)
any(duplicated(ewas_sig_all[, c("CpG")]))
sum(duplicated(ewas_sig_all[, c("CpG")]))

# downloaded the cis-eQTMs files
MA_eQTMs_M1.ga.all <- read.delim("db/eQTMs/MA_eQTMs_M1.ga.all.ann_ov.txt", header = TRUE, stringsAsFactors = FALSE)
colnames(MA_eQTMs_M1.ga.all)
head(MA_eQTMs_M1.ga.all)

eQTMs_all <- left_join(ewas_sig_all, MA_eQTMs_M1.ga.all, by = c("CpG" = "CpG_eQTM"))
write.csv(eQTMs_all, "/results/Downstream/Annotation/eQTMs/eQTMs_all.csv", row.names = FALSE)
# We only used all sample database, considering that subsamples with lower number
# MA_eQTMs_M1.ga.female <- read.delim("db/eQTMs/MA_eQTMs_M1.ga.female.ann_ov.txt", header = TRUE, stringsAsFactors = FALSE)
# MA_eQTMs_M1.ga.male <- read.delim("db/eQTMs/MA_eQTMs_M1.ga.male.ann_ov.txt", header = TRUE, stringsAsFactors = FALSE)
################################################################################


######################
# 8) Significant cis-eQTMs  
######################
eQTMs_all <- read.csv("/results/Downstream/Annotation/eQTMs/eQTMs_all.csv")

#### Significant: pval < 1e-5
eQTMs_sig <- eQTMs_all %>%
  subset(pvalue_eQTM < 1e-5)
#### To convert ENS IDs to Gene symbols
# Connectar a Ensembl
# US East (NCBI)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                host = "https://useast.ensembl.org")

# Llista d'Ensembl IDs
ensembl_ids <- as.character(eQTMs_sig$Gene_eQTM)

# Obtenir gene symbols
result <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = mart
)
eQTMs_final <- eQTMs_sig %>%
  left_join(result, by = c("Gene_eQTM" = "ensembl_gene_id")) 

eQTMs_final <- eQTMs_final[!duplicated(eQTMs_final[, c("CpG", "Gene_eQTM", "Exposure")]), ] 
write.csv(eQTMs_final, file = "/results/Downstream/Annotation/eQTMs/eQTMs_final.csv")
################################################################################
#END

######################
# 9) Annotation for no Gene CpGs  
######################
# library(bumphunter)
# library(minfi)
# library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# library(org.Hs.eg.db)
# library(GenomicRanges)
# Gene location was obtained from Illumina annotation. If the CpG did not have a gene annotated within 1500 bp from the promoter region, the information regarding whether the CpG was downstream or upstream of the gene was obtained using the matchGenes() function from the bumphunter package.
# ewas_sig <- read.csv("/results/Downstream/Annotation/Sig_DMPs_annotated.csv")
# 
# gene_illumina<-ewas_sig$UCSC_RefGene_Name
# 
# #Merge with MatchGenes()
# genome <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
# 
# #Match genes
# ewas_sig$seqnames<-ewas_sig$chr
# ewas_sig$start<-ewas_sig$pos
# ewas_sig$end<-ewas_sig$pos+1
# ewas_sig_gr<-makeGRangesFromDataFrame(ewas_sig, keep.extra.columns=TRUE, start.field="start", end.field="end", seqnames.field="seqnames")
# 
# 
# ewas_sig_genes_bumbhunter<-matchGenes(ewas_sig_gr, genome, type = "fiveprime",  skipExons = TRUE, verbose = TRUE)
# ewas_sig_genes_bumbhunter<-ewas_sig_genes_bumbhunter[,c("name","region","distance","geneL")]
# 
# ewas_sig_final<-cbind(ewas_sig,ewas_sig_genes_bumbhunter)
# # ewas_sig_final<-ewas_sig_final[,c(2:7,9:11,29,31:34,36,37,47:50,54,57,59,85:88)]
# write.table(ewas_sig_final,"Downstream/Annotation/ewas_sig_annotation_genes_illumina_bumbh.txt",sep="\t")

