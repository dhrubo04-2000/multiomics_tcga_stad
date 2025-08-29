setwd("D:/Multiomics_STAD/STAD_multiomics")
library(readr)
library(tidyverse)
STAD_miRNA_minmax <- read_csv("miRNA_data/STAD_miRNA_minmax.csv")
STAD_mRNA_minmax_top_2000 <- read_csv("mRNA_data/STAD_mRNA_minmax_top_2000.csv")
STAD_rppa_minmax <- read_csv("rppa_data/STAD_rppa_minmax.csv")
stad_meth_top2000 <- readRDS("methylation_data/STAD_meth_top2000.rds")

nrow(stad_meth_top2000) 

# Calculate SD across samples
meth_sd <- apply(stad_meth_top2000, 1, sd, na.rm = TRUE)

# Select top 2000 most variable CpGs
top2000_idx <- order(meth_sd, decreasing = TRUE)[1:2000]
stad_meth_top2000 <- stad_meth_top2000[top2000_idx, ]










# Set first column as rownames
STAD_miRNA_minmax <- as.data.frame(STAD_miRNA_minmax)
rownames(STAD_miRNA_minmax) <- STAD_miRNA_minmax[[1]]  # take first column
STAD_miRNA_minmax <- STAD_miRNA_minmax[ , -1]          # drop first column
view(STAD_miRNA_minmax)


#set first column to rownames
STAD_mRNA_minmax_top_2000 <- as.data.frame(STAD_mRNA_minmax_top_2000)
rownames(STAD_mRNA_minmax_top_2000) <- STAD_mRNA_minmax_top_2000[[1]]
STAD_mRNA_minmax_top_2000 <- STAD_mRNA_minmax_top_2000[ , -1]
view(STAD_mRNA_minmax_top_2000)
#set first column to rownames
STAD_rppa_minmax <- as.data.frame(STAD_rppa_minmax)
rownames(STAD_rppa_minmax) <- STAD_rppa_minmax[[1]]
STAD_rppa_minmax <- STAD_rppa_minmax[ , -1]
view(STAD_rppa_minmax)




ncol(STAD_miRNA_minmax)
ncol(STAD_mRNA_minmax_top_2000)
ncol(STAD_rppa_minmax)
ncol(stad_meth_top2000)
# 1. Get the sample IDs from each dataset
colnames(STAD_mRNA_minmax_top_2000)[1:5]
colnames(STAD_miRNA_minmax )[1:5]
colnames(stad_meth_top2000)[1:5]
colnames(STAD_rppa_minmax)[1:5]


# Function to trim TCGA barcodes to first 4 fields: TCGA-XX-YYYY-ZZ
fix_tcga_ids <- function(ids) {
  ids <- gsub("^_", "", ids)  # remove leading underscores
  ids <- gsub("^\\.\\.\\.1$", "ID", ids) # fix weird column "...1" as "ID"
  ids <- sapply(strsplit(ids, "-"), function(x) paste(x[1:4], collapse = "-"))
  return(ids)
}


# Apply to each dataset
colnames(STAD_mRNA_minmax_top_2000) <- fix_tcga_ids(colnames(STAD_mRNA_minmax_top_2000))
colnames(STAD_miRNA_minmax) <- fix_tcga_ids(colnames(STAD_miRNA_minmax))
colnames(stad_meth_top2000) <- fix_tcga_ids(colnames(stad_meth_top2000))
colnames(STAD_rppa_minmax) <- fix_tcga_ids(colnames(STAD_rppa_minmax))



# Now find common samples
common_samples <- Reduce(intersect, list(
  colnames(STAD_mRNA_minmax_top_2000),
  colnames(STAD_miRNA_minmax),
  colnames(stad_meth_top2000),
  colnames(STAD_rppa_minmax)
))

length(common_samples)


mRNA_aligned <- STAD_mRNA_minmax_top_2000[, common_samples]
miRNA_aligned <- STAD_miRNA_minmax[, common_samples]
meth_aligned <- stad_meth_top2000[, common_samples]
RPPA_aligned <- STAD_rppa_minmax[, common_samples]
ncol(mRNA_aligned)




# Optional: sanity check
all(colnames(mRNA_aligned) == colnames(miRNA_aligned))  # should be TRUE
all(colnames(mRNA_aligned) == colnames(meth_aligned))
all(colnames(mRNA_aligned) == colnames(RPPA_aligned))

#save each file as .rds
saveRDS(miRNA_aligned , file = "aligned_samples/STAD_miRNA_minmax_mito.rds")
saveRDS(mRNA_aligned, file = "aligned_samples/STAD_mito_mRNA_minmax_top_2000.rds")
saveRDS(RPPA_aligned, file = "aligned_samples/STAD_mito_rppa_minmax.rds")
saveRDS(meth_aligned, file = "aligned_samples/STAD_meth_mito_top2000.rds")




rownames(mRNA_aligned) <- paste0("mRNA_", rownames(mRNA_aligned))
rownames(miRNA_aligned) <- paste0("miRNA_", rownames(miRNA_aligned))
rownames(meth_aligned) <- paste0("METH_", rownames(meth_aligned))
rownames(RPPA_aligned) <- paste0("RPPA_", rownames(RPPA_aligned))
view(mRNA_aligned)




# After stacking features
multiomics_matrix <- rbind(
  mRNA_aligned,
  miRNA_aligned,
  methylation_aligned,
  RPPA_aligned
)
view(multiomics_matrix)
# Transpose: samples as rows, features as columns
multiomics_matrix_T <- t(multiomics_matrix)

dim(multiomics_matrix_T)  
# Now: n_samples × n_features

#save the final matrix as .csv
write.csv(multiomics_matrix_T, "multiomics_matrix/multiomics_matrix_for_AE.csv", row.names = TRUE)
#save as rds
saveRDS(multiomics_matrix_T, "multiomics_matrix/multiomics_matrix_for_AE.rds")



mutation_maf <- readRDS("mutation_maf.rds")
view(mutation_maf)







