if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")

library(SummarizedExperiment)
library(tidyverse)
library(DESeq2)
library(TCGAbiolinks)

# getting clinical data for TCGA-BRCA cohort -------------------
clinical_brca <- GDCquery_clinic("TCGA-BRCA")
table(clinical_brca$ajcc_pathologic_stage)

# get gene expression data -----------

# build a query to get gene expression data for entire cohort
query_brca_all = GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification",
  sample.type = c("Primary Tumor", "Solid Tissue Normal"),
  access = "open",
  data.format = "TSV")

output_brca <- getResults(query_brca_all)

# download data
GDCdownload(query_brca_all)
# get counts
tcga_brca_data <- GDCprepare(query_brca_all, summarizedExperiment = TRUE)
brca_matrix <- assay(tcga_brca_data, "fpkm_unstrand")

# extract gene and sample metadata from summarizedExperiment object
gene_metadata <- as.data.frame(rowData(tcga_brca_data))
coldata <- as.data.frame(colData(tcga_brca_data))

# remove list columns
coldata_clean <- coldata[, !sapply(coldata, is.list)]

write.csv(brca_matrix, "C:/Users/vanille.lejal/OneDrive - University of Luxembourg/Documents/TT/fastcore_workflow/scr/dataForTesting/brca_matrix.csv", row.names = TRUE)
write.csv(coldata_clean, "C:/Users/vanille.lejal/OneDrive - University of Luxembourg/Documents/TT/fastcore_workflow/scr/dataForTesting/samples_metadata.csv", row.names = TRUE)
write.csv(gene_metadata, "C:/Users/vanille.lejal/OneDrive - University of Luxembourg/Documents/TT/fastcore_workflow/scr/dataForTesting/gene_metadata.csv", row.names = TRUE)


