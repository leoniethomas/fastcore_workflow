library(dplyr)
#library(Seurat)
#library(patchwork)
#library(R.matlab)
#library(tidyverse)
#library(DESeq2)
#library(devtools)
library(biomaRt) #BiocManager::install("biomaRt")
require(gridExtra)
library(grid)
library(Biobase)
#library(sva)
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


counts_to_tpm <- function(count_matrix, gene_lengths) {
  # Convert gene lengths to kilobases
  gene_lengths_kb <- gene_lengths / 1000
  
  # Calculate RPK for each sample
  rpk_matrix <- sweep(count_matrix, 1, gene_lengths_kb, "/")
  
  # Calculate scaling factor for each sample
  scaling_factors <- colSums(rpk_matrix) / 1e6
  
  # Calculate TPM for each sample
  tpm_matrix <- sweep(rpk_matrix, 2, scaling_factors, "/")
  
  return(tpm_matrix)
}

#FPKM <- apply(count_df, 2, function(x) counts_to_fpkm(x, lengths ))

get_gene_length <- function(count_df, identifier_gene_in_data = "hgnc_symbol" ) { #  change in the code 
  # ensembl_gene_id: based on the input gene format ensembl_gene_id/ hgnc_symbol
  # Name: colname of the gene ids
  df <- count_df
  geneSymbols <-  df$Name

  human <- useMart("ensembl", dataset="hsapiens_gene_ensembl") # hgnc_symbol
  gene_coords=getBM(attributes=c(identifier_gene_in_data,"ensembl_gene_id", "start_position","end_position"), 
                    filters=identifier_gene_in_data, values=geneSymbols, mart=human)
  colnames(gene_coords)[colnames(gene_coords) == identifier_gene_in_data] <- "Name"
  gene_coords$size=gene_coords$end_position - gene_coords$start_position
  # Taking the maximum length for dublicated lengths
  gene_coords %>% group_by(Name) %>% summarise(size= max(size)) -> gene_coords
  gene_coords <- gene_coords[gene_coords$Name %in% df$Name,]
  # Taking the maximim of count for dublicated genes
  df %>% group_by(Name) %>% summarize(across(everything(), list(max)) )-> df
  
  df_knw_len <- df[df$Name %in% gene_coords$Name,]
  
  df_len <- left_join(gene_coords,df_knw_len, by = "Name")
  df_len <- df_len[!is.na(df_len$size),]
  return(df_len)
}

### parameters
# is this data already normalized, these are no raw counts

wd_path = "/Users/leonie.thomas/20250225_glynn_bulk_metabolic_model/"
output_path = paste0(wd_path, "data/bulkRNAseq/")

file_name = "data/bulkRNAseq/Merged_raw_counts.txt"
#file_name_metadata = "data/GSE176078_Wu_etal_2021_bulkRNAseq_metadata"

### load counts and get gene length

counts <- read.table(	paste0(wd_path, file_name), sep = "\t", #skip = 1,
			header = TRUE, row.names = 1)

counts$Name <- rownames(counts)
counts <- get_gene_length(counts, identifier_gene_in_data = "hgnc_symbol" )
write.table(counts, paste0(output_path , "raw_counts_with_genelengths.csv"), quote = FALSE, row.names = FALSE)
## get the fpkm 

df_fpkm <- as.data.frame(apply(counts[,3:ncol(counts)], 2, function(x) counts_to_fpkm(x, counts$size)))
rownames(df_fpkm) <- counts$Name

colnames(df_fpkm) <- gsub("_1$", "", colnames(df_fpkm))
write.table(df_fpkm, paste0(output_path , "data_fpkm.csv"), quote = FALSE)



### calculate tpm 

#df_tpm <- counts_to_tpm(counts[,3:ncol(counts)],counts$size)

# %colSums(df_tpm)
# %colnames(df_tpm) <- gsub("_1$", "", colnames(df_tpm))
# %write.table(df_tpm, paste0(output_path , "data_tpm.csv"), quote = FALSE)
# 
# 
# ### perform pca on the samples
# 
# 
# %results <- prcomp(t(df_fpkm))
# 
# %results$rotation <- -1*results$rotation
# 
# %results$x %>%
# %%  ggplot( aes(x=PC1, y=PC2)) +
# %%  geom_point(size = 4) +
# %%  labs(x = paste0("PC1 (expl. Var: ",as.character(summary(results)$importance[2,1])," )"),
# %%       y = paste0("PC2 (expl. Var: ",as.character(summary(results)$importance[2,2])," )"))
# 
# 
