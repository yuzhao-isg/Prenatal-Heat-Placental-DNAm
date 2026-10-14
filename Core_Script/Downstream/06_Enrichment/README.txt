# Pathway Enrichment Analysis

## 1. Overview

This folder contains the enrichment analysis scripts for the overall sample and sex-stratified analyses of the epigenome-wide association study (EWAS).

The analytical workflow for the sex-stratified analyses is the same as that used for the overall-sample analyses. The only difference is the input file path, which points to the quality-controlled (QC) results for the corresponding sex-stratified analysis. The same analysis procedures and settings are then applied to the sex-stratified results.

## 2. Overall-Sample and Sex-Stratified Analyses

The sex-stratified enrichment analyses follow the same workflow as the overall-sample analyses:

1. Use the QC-passed EWAS results for the relevant analysis.
2. Update the input file path to point to the corresponding sex-stratified results.
3. Run the same enrichment analysis scripts using the same analytical procedures and settings.

No other changes to the analysis workflow are required, provided that the input files have the expected format and contain the required variables.

## 3. Enrichment Analysis Methods

Two complementary approaches are used to investigate biological pathways and functional categories associated with the EWAS findings: Gene Set Enrichment Analysis (GSEA) and Over-Representation Analysis (ORA).

### 3.1. GSEA

GSEA is a rank-based enrichment approach that evaluates whether genes belonging to a predefined gene set are concentrated toward either end of a ranked gene list.

Rather than selecting genes using a predefined significance threshold, GSEA uses information from the broader set of results to assess enrichment. The ranking metric and gene-mapping procedure depend on the specific implementation used in the analysis.

### 3.2. ORA

ORA evaluates whether a set of selected genes is over-represented in particular biological pathways or functional categories compared with a reference gene background.

In this analysis, CpGs are selected using a predefined P-value threshold (P < 1 × 10⁻⁵). The selected CpGs are then mapped to genes, which are used as input for the enrichment analysis.

Unlike GSEA, ORA relies on a threshold-defined set of genes rather than the full ranked gene list.

## 4. Key Differences Between GSEA and ORA

| Feature             | GSEA                                                             | ORA                                                                              |
| ------------------- | ---------------------------------------------------------------- | -------------------------------------------------------------------------------- |
| Input               | Ranked gene list                                                 | Selected genes derived from threshold-filtered results                           |
| Selection threshold | No predefined significance threshold for selecting genes         | P < 1 × 10⁻4 for selecting CpGs in this analysis                                 |
| Main objective      | Identify gene sets enriched toward either end of the ranked list | Identify pathways or functional categories over-represented among selected genes |
| Use of information  | Uses the broader ranked list                                     | Focuses on genes passing the selection threshold                                 |
| Interpretation      | Captures coordinated patterns across the ranked results          | Identifies enrichment among the selected findings                                |

## 5. Notes

- The overall-sample and sex-stratified analyses use the same analytical workflow.
- For sex-stratified analyses, ensure that the input paths point to the appropriate QC-passed results for each sex.
- GSEA and ORA are complementary approaches and may identify different biological pathways because they use different input strategies.
- The P-value threshold used for ORA is a selection criterion for downstream enrichment analysis and should not be interpreted as a genome-wide significant threshold for the EWAS itself.
- Gene annotation, CpG-to-gene mapping, and the background gene set should follow the specifications of the corresponding enrichment method and software implementation.
