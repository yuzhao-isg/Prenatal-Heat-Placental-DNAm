#-----------------------------------------------------------------------------#
#              Enrichment analysis - fetal placenta chromatin states          #
# Based on ROADMAP data:                        
#   - FP_15_E091: 15 chromatin states
#   - FP_18_E091: 18 chromatin states
#-----------------------------------------------------------------------------#
rm(list = ls())
setwd("~/010_EWAS_temp_YZ") # change to your path

######################
# 1) Load packages
######################
library(EASIER)
library(GenomicRanges)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(dplyr)
library(readr)
library(stringr)
library(purrr)

anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

######################
# 2) Prepare CpG files
######################
base_dir <- "/results/Main_Model/CpG_quality_control/"
files <- list.files(
  base_dir,
  pattern = "Exclude_CpGs_results_.*\\.RData$",
  full.names = TRUE
)

# Get the exposure name
get_exposure <- function(path) {
  fname <- basename(path)
  fname <- sub("Exclude_CpGs_results_", "", fname)
  fname <- sub("\\.RData$", "", fname)
  return(fname)
}

results <- list()


######################
# 3) Function
######################
#' Find Overlaps between two genomics ranges and returns the coincident values. Used to search overlaps between CpGs (selected data) and Fetal Placenta (complete data)
#'
#' @param selected Genomic Ranges Variable that we are interested to enrich with data obtained from 'comlete' genomic ranges
#' @param complete Genomic Ranges with Fetal Placenta data 15 or 18 or other genomic ranges that overlaps with selected genomic ranges
#'
#' @return A list with the Overlap results
#' \itemize{
#'  \item{"ranges"}{A subset by overlaps}
#'  \item{"hits"}{Coincident data between selected and complete}
#'  \item{"shits"}{Subject hits}
#'  \item{"qhits"}{Query hits}
#'  \item{"values"}{Coincident names associated to each position}
#' }
#'
#' @export
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
######################  
#' Fisher test
#'
#' Get all Fisher tests taking in to account a significative variable parameter like FDR values. This function gets all chromatin state regressions and optionally store results in a file
#'
#' @param significative numerical. Variable to take in to account as significative variable, could be FDR, p-value,...
#' @param position dataframe with all relative positions to island values to perform regression
#' @param outputdir optional string. Output path to store file with results, by default results are written in current dir
#' @param outputfile optional string. File name to store results if no name is provided results are not written.
#' @param plots boolean. If plot is TRUE, plot results
#'
#' @return
#'
#' @export
getAllFisherTest <- function(significative, position, outputdir = ".", outputfile = NULL, plots = TRUE) {
  positions <- getUniqueValues(position)

  lregs <- lapply(
    positions,
    function(x) getFisherTest(significative, ifelse(position == x, "yes", "no"), x)
  )

  ans <- data.frame(matrix(unlist(lregs), nrow = length(positions), byrow = TRUE))
  colnames(ans) <- names(lregs[[1]])

  if (!is.null(outputfile)) {
    if (!is.null(outputdir) && !is.na(outputdir) && outputdir != ".") {
      dir.create(outputdir, showWarnings = FALSE, recursive = TRUE)
    }

    filename <- tools::file_path_sans_ext(basename(outputfile))
    write.csv(ans, paste0(file.path(outputdir), "/Fisher_", filename, ".csv"))
  }

  if (plots) {
    plot_OR(ans, outputdir = outputdir, outputfile = outputfile)
  }

  return(ans)
}
##############################
#' Get Fisher test for significative variable attending to criteria
#'
#' Get Fisher test by criteria taking in to account a significative variable parameter like FDR values
#'
#' @param significative numerical. Variable to take in to account as significative variable, could be FDR, p-value,...
#' @param criteria vector with criteria
#' @param varname string. criteria Relative to Island name. For example : "Island", "OpenSea", "S_Shore", "N_Shore", "S_Shelf", "N_Shelf" or other criteria like gene position
#'
#' @return
#'
#' @export
getFisherTest <- function(significative, criteria, varname)
{
  
  rp <- as.data.frame.matrix(table(significative, criteria))
  
  if( all(rownames(rp) %in% c('yes', 'no')) ){
    if(dim(rp)[1]<2) {
      rp <- rbind(rp,c(0,0))
      if(!'yes' %in% rownames(rp)) {
        rownames(rp) <- c('no', 'yes')
      }else {
        rownames(rp) <- c('yes', 'no')
      }
    }
  }
  
  if( any(rownames(rp) %in% c('Hypo-yes','Hyper-yes' )) &&
      sum(rownames(rp) %in% c('Hypo-yes','Hyper-yes','Hypo-no','Hyper-no' ))<4 ) {
    rnames <- rownames(rp)
    rp <- rbind(rp,c(0,0))
    if(!'Hyper-yes' %in% rownames(rp)) {
      rownames(rp) <- c(rnames, 'Hyper-yes')
    }
    if(!'Hypo-yes' %in% rownames(rp)) {
      rownames(rp) <- c(rnames, 'Hypo-yes')
    }
    if(!'Hyper-no' %in% rownames(rp)) {
      rownames(rp) <- c(rnames, 'Hyper-no')
    }
    if(!'Hypo-no' %in% rownames(rp)) {
      rownames(rp) <- c(rnames, 'Hypo-no')
    }
  }
  
  
  if(dim(rp)[2]<2) {
    rp <- cbind(rp,c(0,0))
    if(!'yes' %in% colnames(rp)) {
      rp$yes <- rep(0,2)
    }else {
      rp$no <- rep(0,2)
    }
  }
  
  rp <- rp[,c("yes","no")]
  
  if(all(rownames(rp) %in% c('yes', 'no'))){
    rp1 <- rbind(nosig = rp[1, ], sig = colSums(rp[-1, ]))[c(2,1),]
    xt <- chisq.test(rp1)
    or <- (rp1[1, 1] * rp1[2, 2]) / (rp1[1, 2] * rp1[2, 1])
  } else {
    rp1 <- rp[c("Hyper-yes","Hypo-yes"), ]
    xt <- chisq.test(rp1)
    or <- (rp1[1, 1] * rp1[2, 2]) / (rp1[1, 2] * rp1[2, 1])
  }
  
  return(list("Data" = varname,
              "OR" = or,
              "OR.inf" = exp(log(or)-1.96*sqrt(sum(1/rp1))),
              "OR.sup" = exp(log(or)+1.96*sqrt(sum(1/rp1))),
              "p-val" = xt$p.value))
}
  
################################################################################


######################
# 4) Run enrichment for each exposure
######################
for (f in files) {
  message("Processing: ", f)
  
  # Read data
  obj_name <- load(f)
  df <- get(obj_name)
  
  # Exposure name
  exp <- get_exposure(f)
  anno_df <- as.data.frame(anno)
  data <- merge(
    df,
    anno_df[, c("Name", "chr", "pos")],
    by.x = "CpG",
    by.y = "Name",
    all.x = TRUE
  )
  
  # Convert to GenomicRanges
  data.GRange <- GRanges(
    seqnames = Rle(data$chr),
    ranges = IRanges(data$pos, end = data$pos),
    name = data$CpG,
    chr = data$chr,
    pos = data$pos
  )
  names(data.GRange) <- data.GRange$name
  
  # Overlap with fetal placenta chromatin states
  over15 <- findOverlapValues(data.GRange, FP_15_E091)
  
  # Add chromatin-state annotations to CpGs
  data.chrstates <- c(mcols(over15$ranges), over15$values)
  colnames(data.chrstates)[grep("States", colnames(data.chrstates))] <-
    c("States15_FP")
  
  # Merge annotated data with chromatin states
  crom_data <- merge(data, data.chrstates, by.x = "CpG", by.y = "name")
  
  # Hypo / hyper split
  crom_data$signif_hypo <- ifelse(crom_data$pval < 1e-4 & crom_data$coef < 0, "yes", "no")
  crom_data$signif_hyper <- ifelse(crom_data$pval < 1e-4 & crom_data$coef > 0, "yes", "no")
  
  # Fisher test for hyper- and hypo-methylated CpGs
  States15FP_hyper <- getAllFisherTest(
    crom_data$signif_hyper,
    crom_data$States15_FP,
    outputdir = "/results/Downstream/Enrichment/CromStates/Fisher_Hyper",
    outputfile = exp
  )
  
  States15FP_hypo <- getAllFisherTest(
    crom_data$signif_hypo,
    crom_data$States15_FP,
    outputdir = "/results/Downstream/Enrichment/CromStates/Fisher_Hypo",
    outputfile = exp
  )
  
  # Plot collapsed results
  plot_TestResults_Collapsed(
    list(hypo = States15FP_hypo, hyper = States15FP_hyper),
    outputdir = "/results/Downstream/Enrichment/CromStates",
    outputfile = exp
  )
}
#######################################


######################
# 5) Combine Fisher test results
######################
root_dir <- "/results/Downstream/Enrichment/CromStates"

csv_files <- list.files(
  path = root_dir,
  pattern = "\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

all_data <- map_dfr(csv_files, function(file) {
  expression_level <- basename(dirname(file))
  fname <- tools::file_path_sans_ext(basename(file))
  tmp <- str_remove(fname, "^Fisher_")

  group <- case_when(
    str_detect(tmp, "_F$") ~ "Female",
    str_detect(tmp, "_M$") ~ "Male",
    TRUE ~ "overall"
  )

  exposure <- tmp %>%
    str_remove("_F$") %>%
    str_remove("_M$")

  dat <- read_csv(file, show_col_types = FALSE)

  dat <- dat %>%
    mutate(FDR_p = p.adjust(`p-val`, method = "BH"))

  dat %>%
    mutate(
      Expression_level = expression_level,
      exposure = exposure,
      group = group
    )
})

all_data <- all_data %>%
  mutate(
    exposure = recode(
      exposure,
      "indoor_HI95" = "Indoor Heat Index",
      "indoor_T95" = "Indoor Temperature",
      "outdoor_HI95" = "Outdoor Heat Index",
      "outdoor_T95" = "Outdoor Temperature"
    )
  )

write_csv(all_data, file.path(root_dir, "CromStates_enrich_results.csv"))
################################################################################
# END