#======================================================================
# Update 13.10.2026 by Yu Zhao
# original prepared by Sofia Aguilar Lacasana
# Final - RUN EWAS PRENATAL HEAT EXPOSURE VS PLACENTAL DNA METHYLATION
#======================================================================
rm(list = ls())
######################
# 1) Load packages
######################
library(PACEanalysis)
library(data.table)
library(parallel)
library(dplyr)

######################################################
# 1) Load phenodata with the exposures and covariates
######################################################
setwd("~/010_EWAS_temp_YZ")
load("/db/phenodataframe.rds")
phenodataframe <- phenodataframe %>%
  subset(trimester == "Whole_pregnancy")

####################################################################
#2) Load methylation data processed data and cell type proportions
####################################################################

#---------------------------------------------------------
#Load the beta matrix after the QC and without outliers
#----------------------------------------------------------
load("/db/BISC_Pla_MethylEPIC_PreprocessedBetas_nooutliers_20231115.RData")
#Keep the same number of samples than in the pheno
class(betafinal.nooutlier)
betafinal.nooutlier<-as.data.frame(betafinal.nooutlier)
head(colnames(betafinal.nooutlier))
samples.included <- as.character(unlist(phenodataframe[, 3]))
betafinal.nooutlier.inc <- betafinal.nooutlier[,colnames(betafinal.nooutlier) %in% samples.included]
dim(betafinal.nooutlier.inc)
#---------------------------------------------------------------------
#Load these preprocessed object which has the celltype proportions
#---------------------------------------------------------------------
load("/db/BISC_Pla_MethylEPIC_Preprocessed_20231115.RData")
celltypes<-processedOut$Omega
class(celltypes)
head(celltypes)
celltypes<-as.data.frame(celltypes)
head(rownames(celltypes))
#Keep the same number of samples than in the pheno
celltypes.inc<-celltypes[rownames(celltypes) %in% samples.included,]
dim(celltypes.inc)


###########################################################
#3) Check order across the samples of the different objects 
###########################################################
#--------------------------
#pheno + matrix DNA meth
#--------------------------
table(ifelse(phenodataframe$Basename==colnames(betafinal.nooutlier.inc),"Matched","--NOT MATCHED--"))
# Order samples
phenodataframe1<-phenodataframe[order(match(phenodataframe$Basename,colnames(betafinal.nooutlier.inc))),]
table(ifelse(phenodataframe1$Basename==colnames(betafinal.nooutlier.inc),"Matched","--NOT MATCHED--"))

#-----------------
#pheno + celltypes
#-----------------
table(ifelse(phenodataframe1$Basename==rownames(celltypes.inc),"Matched","--NOT MATCHED--"))

celltypes.inc<-as.matrix(celltypes.inc)
dim(celltypes.inc)
head(celltypes.inc)

betafinal.nooutlier.inc<-as.matrix(betafinal.nooutlier.inc)
dim(betafinal.nooutlier.inc)
head(betafinal.nooutlier.inc)

#############################
#5 Run analysis
#############################

## if running in parallel, checking the number of available cores 

#-----------------
#-----------------
# SENSITIVITY ANALYSIS - only run for Europeans
#-----------------
#-----------------
allvarsofinterest=c("count_HI_exposure_95_indoor_day","count_heat_exposure_95_indoor_day",
                    "count_HI_exposure_95_home","count_heat_exposure_95_home")

colnames(phenodataframe1)[colnames(phenodataframe1) %in% c("sex_0y_c")] <- "Sex"
phenodataframe1$Sex[phenodataframe1$Sex == "Boy"] <- "Male"
phenodataframe1$Sex[phenodataframe1$Sex == "Girl"] <- "Female"

for (i in 1:length(allvarsofinterest)){
  
  cat("Exposure:",allvarsofinterest[i],"\n")
  tempresults.MAIN<-dataAnalysis(phenofinal=phenodataframe1, #add your final phenodataframe
                                 betafinal=betafinal.nooutlier.inc, #add your final methylation matrix. DEMO first with 100 CpGs --> CHNAGE LATER AND RE-RUN
                                 array="EPIC", #specify the array type
                                 maxit=100, #default
                                 robust=TRUE, # specify "TRUE" if you want to run robust linear regression models
                                 Omega=celltypes.inc, #add your final celltype db
                                 vartype="ExposureCont", #specify which type are the variables of interest - in this case, we are interested in continuous exposures
                                 varofinterest=allvarsofinterest[i],
                                 #Add all variables that you want to add in the descriptive table
                                 Table1vars=c("count_HI_exposure_95_indoor_day","count_heat_exposure_95_indoor_day","count_HI_exposure_95_home",
                                              "count_heat_exposure_95_home","Sex","age_12w_m","educ_level_m_2cat",
                                              "ses_income_acu","hospital_delivery_3c", "Meanlog2oddsContamination","Active_smok_any",
                                              "gestage_0y_c_weeks","parity_m_2cat","conception_month","covid_frontier_m_3c"),
                                 StratifyTable1=FALSE, #in case you want to create de descriptive table stratifying by an specific variable
                                 StratifyTable1var=NULL,
                                 #Add the adjusment variables (covariates)
                                 adjustmentvariables=c("Sex","age_12w_m","educ_level_m_2cat","ses_income_acu","hospital_delivery_3c",
                                                       "Meanlog2oddsContamination","Active_smok_any","gestage_0y_c_weeks","parity_m_2cat","ns_cd_1","ns_cd_2", 
                                                       "ns_cd_3","ns_cd_4","ns_cd_5","ns_cd_6","ns_cd_7","ns_cd_8","ns_cd_9","ns_cd_10","ns_cd_11","ns_cd_12"),# "conception_month","covid_frontier_m_3c"
                                 RunUnadjusted=FALSE,#model without adjusting for any covariate
                                 RunAdjusted=TRUE,#model without adjusting for celltype proportions
                                 RunCellTypeAdjusted=TRUE,#model adjusted for celltype proportions (in my case, the MAIN model)
                                 RunSexSpecific=TRUE,# run analyses stratified by sex
                                 RunCellTypeInteract=TRUE,#run interaction by celltype proportion
                                 RestrictToSubset=TRUE,
                                 RestrictionVar="ethnicity_c_2cat",
                                 RestrictToIndicator="eur",
                                 RestrictionVar=NULL,
                                 RestrictToIndicator=NULL,
                                 number_cores=4, #default
                                 runparallel=TRUE,#default
                                 destinationfolder="/results/Sensitivity_Model/",
                                 savelog=TRUE,#default
                                 cohort="BISC",analysisdate="20261013",
                                 analysisname="MAIN")
  
}
