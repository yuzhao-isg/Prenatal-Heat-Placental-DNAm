# Prenatal heat exposure and placental DNA methylation: EWAS in the BiSC cohort

This repository contains the analysis code for an epigenome-wide association study (EWAS) of prenatal indoor and outdoor heat exposure and placental DNA methylation in the [Barcelona Life Study Cohort (BiSC)](https://projectebisc.org/en/home/).

The workflow covers preparation of the analytical dataset, robust linear regression in the overall sample and sex-stratified samples, CpG probe filtering, gene and regulatory annotation, enrichment analyses, comparison with previous EWAS findings, and generation of study figures.

## Study overview

| Item | Description |
|---|---|
| Study population | 565 placental samples from the BiSC cohort, Barcelona, Spain |
| DNA methylation platform | Illumina Infinium MethylationEPIC array |
| Outcome | Placental DNA methylation beta-values, ranging from 0 to 1 |
| Exposure | Cumulative prenatal heat-exposure days, assessed separately for indoor and outdoor temperature and heat index (HI) |
| Main exposure threshold | Calendar-day-specific 95th percentile |
| Statistical model | Robust linear regression |
| Analysis populations | Overall sample, females, and males |

## Repository organisation

The repository is organised into four main directories:

| Directory | Contents |
|---|---|
| `Core_Script/` | Core EWAS scripts and the `Downstream/` subdirectory for annotation and follow-up analyses |
| `db/` | Reference data for imprinting-region annotation |
| `Package_Docs/` | Documentation and usage guidelines for R packages used in the analyses |
| `Plot/` | Scripts for generating the study figures |

### Core_Script

| File | Purpose |
|---|---|
| `01_Prepare_database.R` | Prepares the analytical dataset containing heat-exposure measures, placental DNA methylation outcomes, and covariates. |
| `02_EWAS_Analysis` | Runs the main EWAS using robust linear regression in the overall sample and in female and male strata. |
| `03_CpGs_QC.R` | Applies the specified probe exclusions to define the CpG set retained for reporting and downstream analyses. |

### Core_Script/Downstream

| File or directory | Purpose |
|---|---|
| `04_Gene_Annotation.R` | Maps CpGs to genes using Illumina positional annotations, trophoblast enhancer–gene links from EPIraction, and placental cis-eQTMs from the PACE catalog. |
| `05_Annotation_Regulatory_Elements.R` | Annotates CpGs using placental cis-mQTLs, partially methylated domains (PMDs), germline differentially methylated regions (gDMRs), chromatin states, and imprinting control regions (ICRs). |
| `06_Enrichment/` | Contains scripts for functional enrichment and enrichment in placental genomic features; individual scripts are described below. |
| `07_EWAS_Overlap.R` | Queries EWAS Atlas and the EWAS Catalog to identify previously reported associations involving significant CpGs from this study. |
| `08_Tissuecell_gene_protein_expression.txt` | Provides the workflow for characterising genes mapped to significant CpGs or linked through significant eQTM/enhancer findings, including tissue and cell-type expression, protein expression, placental enrichment, and potential transcriptional regulation in the placenta. |

#### Enrichment scripts

Location: `Core_Script/Downstream/06_Enrichment/`

| File | Analysis |
|---|---|
| `GSEA_Overall.R` | Performs functional enrichment analyses using the `MSigDB_Hallmark_2020`, `DisGeNET`, and `ENCODE_and_ChEA_Consensus_TFs` gene-set libraries, covering biological pathways, disease associations, and transcription-factor target sets. |
| `CromStates_Placenta.R` | Evaluates enrichment of CpGs with EWAS P values below `P < 1e-4` in fetal placental chromatin states. |
| `PMD_Placenta.R` | Evaluates enrichment of CpGs with EWAS P values below `P < 1e-4` in placental PMDs, using the regions reported by Schroeder et al. (2013). |

### db

Contains the `imprinting_regions` reference resource used for imprinting-related annotation.

Other external reference resources are specified in the scripts that use them. Their file paths must match the local analysis environment.

### Package_Docs

Contains documentation and usage guidelines for R packages used in the workflow. These documents support the implementation and interpretation of the relevant analysis steps.

### Plot

| File | Figure or comparison |
|---|---|
| `Enrichment_Heatmap.R` | Generates heatmaps of enrichment results for placental chromatin states and placental PMDs. |
| `Manhattan_Volcano Plot.R` | Generates Manhattan and volcano plots of the EWAS results. |
| `QC_Plot.R` | Generates quantile–quantile (QQ) plots and P-value distribution plots. |
| `Sensitivity_Main_Correlation_Plot.R` | Compares effect estimates between the main and sensitivity analyses for Bonferroni-significant and suggestive CpGs, visualising the consistency of the associations. |

## Data availability

The companion EWAS summary-statistics dataset contains CpG-level association estimates and genomic coordinates. The exported regression coefficients and standard errors are expressed per 1 additional heat-exposure days on the original methylation beta-value scale.

**Zenodo summary-statistics record:** 10.5281/zenodo.23368399.

## Associated publication and citation

The associated research article has not yet been published. The full citation and article DOI will be added when available.

When using this code, please cite the repository and the associated article once published. If using the deposited summary statistics, also cite the corresponding version of the Zenodo dataset.

