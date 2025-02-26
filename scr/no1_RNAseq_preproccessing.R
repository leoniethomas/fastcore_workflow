library(dplyr)
library(Seurat)
library(patchwork)
library(R.matlab)
library(tidyverse)
library(DESeq2)
library(devtools)
library(biomaRt) #BiocManager::install("biomaRt")
require(gridExtra)
library(grid)
library(Biobase)
library(sva)
#library(ggpubr)
library(tidyverse)
library(tidyr)
require(dplyr)
library(rrcov)
library(data.table)
library(vroom)
library(tibble)
library(grid)
library(vroom)
library(data.table)


# this script was written on the basis of : https://github.com/sysbiolux/ISB705MetabolicNetworkModeling/blob/main/DataProcessing/R/Count_to_FPKM_Example.R


counts_to_fpkm <- function(counts, lengths) {
  exp(log(counts) + log(1e9) - log(lengths) - log(sum(counts)) )
}

#FPKM <- apply(count_df, 2, function(x) counts_to_fpkm(x, lengths ))

calculate_FPKM<- function(count_df) { #  change in the code 
  # ensembl_gene_id: based on the input gene format ensembl_gene_id/ hgnc_symbol
  # Name: colname of the gene ids
  df <- count_df
  geneSymbols <-  df$Name

  human <- useMart("ensembl", dataset="hsapiens_gene_ensembl") # hgnc_symbol
  gene_coords=getBM(attributes=c("hgnc_symbol","ensembl_gene_id", "start_position","end_position"), 
                    filters="hgnc_symbol", values=geneSymbols, mart=human)
  gene_coords$size=gene_coords$end_position - gene_coords$start_position
  # Taking the maximum length for dublicated lengths
  gene_coords %>% group_by(hgnc_symbol) %>% summarise(size= max(size)) -> gene_coords
  gene_coords <- gene_coords[gene_coords$hgnc_symbol %in% df$Name,]
  # Taking the maximim of count for dublicated genes
  df %>% group_by(Name) %>% summarize(across(everything(), list(max)) )-> df
  
  df_knw_len <- df[df$Name %in% gene_coords$hgnc_symbol,]
  df_len <- left_join(gene_coords,df_knw_len,by=c('hgnc_symbol'="Name"))
  df_len <- df_len[!is.na(df_len$size),]

  df_fpkm <- apply(df_len[,3:ncol(df_len)], 2, function(x) counts_to_fpkm(x, df_len$size))
  df_fpkm <- data.frame(df_fpkm)
  rownames(df_fpkm) <- df_len$hgnc_symbol
  df_fpkm
}


### parameters
# is this data already normalized, these are no raw counts

wd_path = "/Users/leonie.thomas/20241104_TNBC_wu_2021/"

file_name = "data/GSE176078_Wu_etal_2021_bulkRNAseq_raw_counts"
file_name_metadata = "data/GSE176078_Wu_etal_2021_bulkRNAseq_metadata"

### investigate the counts

counts <- read.table(	paste0(wd_path, file_name, ".txt"), sep = "\t", #skip = 1,
			header = TRUE, row.names = 1)

counts$Name <- rownames(counts)
df_fpkm <- calculate_FPKM(counts)
colnames(df_fpkm) <- gsub("_1$", "", colnames(df_fpkm))
write.table(df_fpkm, paste0("./", file_name, "_fpkm.csv"), quote = FALSE)


####

#meta_data <- read.table(paste0(wd_path, file_name_metadata, ".txt"), sep = "\t",
#			header = TRUE, row.names = 1) %>%
#			tibble::rownames_to_column("CaseID") %>%
#			dplyr::mutate(CaseID = gsub("-", "",paste0("CID",CaseID)))
#
#common_caseids = intersect(meta_data$CaseID,colnames(counts))			


#counts <- counts[,common_caseids]
#rownames(meta_data) <- meta_data$CaseID
#meta_data <- meta_data[common_caseids,]

#counts <- counts[rowSums(is.na(counts)) == 0,]


#####
# hugues code

#library(DESeq2)

# Create the DESeq2 object
#dds <- DESeqDataSetFromMatrix(countData = round(counts, digits = 0), 
#                              colData   = meta_data,
#			      design = ~ Subtype.by.IHC) 

#dds <- estimateSizeFactors(dds)
#dds <- dds[rowSums(counts(dds, normalized = FALSE)) > 100, ]
#dds <- DESeq(dds, parallel = TRUE)

# Normalized counts
#vst <- vst(dds, blind = FALSE)


#write.table(counts(dds, normalized = TRUE), paste0("./", file_name, "_normalized_DeSeq.csv"), quote = FALSE)


###

#tnbc_bulk <- CreateSeuratObject(counts = counts,
#				project = "tnbc_bulk",
#				min.cells = 2)

#tnbc_bulk <- NormalizeData(tnbc_bulk,
#				normalization.method = "LogNormalize",
#				scale.factor = 10000)

#pbmc <- ScaleData(pbmc, features = all.genes) #??? do this too ? 

#mat <- GetAssayData(object = tnbc_bulk, assay = "RNA", slot = "data")

#write.table(mat, paste0("./", file_name, "_normalized.csv"), quote = FALSE)
#writeMat(paste0("./", file_name, "_normalized.mat"), labpcexport = mat)


### counts to FPKM

#library(countToFPKM)

#file.readcounts <- system.file("extdata", file_name, package="countToFPKM")
#file.annotations <- system.file("extdata", "Biomart.annotations.hg38.txt", package="countToFPKM")
#file.sample.metrics <- system.file("extdata", "RNA-seq.samples.metrics.txt", package="countToFPKM")

# Import the read count matrix data into R.
#counts <- as.matrix(read.csv(file.readcounts))


# Import feature annotations. 
# Assign feature lenght into a numeric vector.
#gene.annotations <- read.table(file.annotations, sep="\t", header=TRUE)
#featureLength <- gene.annotations$length

# Import sample metrics. 
# Assign mean fragment length into a numeric vector.
#samples.metrics <- read.table(file.sample.metrics, sep="\t", header=TRUE)
#meanFragmentLength <- samples.metrics$meanFragmentLength

# Return FPKM into a numeric matrix.
#fpkm_matrix <- fpkm(counts, featureLength, meanFragmentLength)