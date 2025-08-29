setwd("D:/Multiomics_STAD/STAD_multiomics")
library(readr)
library(tidyverse)
library(readr)
STAD_matrix_mapped <- read_csv("mRNA_data/STAD_matrix_mapped.csv")
View(STAD_matrix_mapped)
nrow(STAD_matrix_mapped)  # Check number of genes
# Remove genes with zero expression in >20% of samples
keep <- rowSums(STAD_matrix_mapped == 0) / ncol(STAD_matrix_mapped) <= 0.2
STAD_mRNA_filtered <- STAD_matrix_mapped[keep, ]
nrow(STAD_mRNA_filtered)  # Check number of genes after filtering)

view(STAD_mRNA_filtered)

#set first column to gene name
colnames(STAD_mRNA_filtered)[1] <- "Genes"
# Set rownames and convert all other columns to numeric in one step
STAD_mRNA_filtered <- as.data.frame(
  sapply(STAD_mRNA_filtered[ , -1], function(x) as.numeric(as.character(x))),
  row.names = STAD_mRNA_filtered$Genes
)
view(STAD_mRNA_filtered)

# Convert to matrix if it's a data frame
STAD_mat <- as.matrix(STAD_mRNA_filtered)
view(STAD_mat)


# Create a dummy condition table (DESeq2 requires colData)
# Replace with your actual metadata if available
sample_info <- data.frame(
  row.names = colnames(STAD_mat),
  condition = rep("A", ncol(STAD_mat)) # placeholder
)
view(sample_info)
# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = STAD_mat,
  colData = sample_info,
  design = ~ 1
)





# Apply VST
vsd <- vst(dds, blind = TRUE)

# Get transformed data
STAD_vst <- assay(vsd)

# View transformed matrix
View(STAD_vst)


# Compute standard deviation for each gene (row)
gene_sd <- apply(STAD_vst, 1, sd)

# Sort genes by decreasing SD
gene_sd_sorted <- sort(gene_sd, decreasing = TRUE)

# Optional: select top N most variable genes, e.g., top 2000
top_genes <- names(gene_sd_sorted)[1:2000]

# Subset the matrix to top variable genes
STAD_vst_top <- STAD_vst[top_genes, ]
view(STAD_vst_top)

#save the file
write.csv(STAD_vst_top, "mRNA_data/STAD_vst_top_2000.csv", row.names = TRUE)

# Load the data(mitochondrial genes)
library(readxl)
Human_MitoCarta3_0 <- read_excel("Human.MitoCarta3.0.xls", 
                                 sheet = "A Human MitoCarta3.0")


#filter for mitochondrial genes
mito_genes <- Human_MitoCarta3_0$Symbol

# Filter rows where rownames match mito_genes
STAD_mito_VST <- STAD_vst[rownames(STAD_vst) %in% mito_genes, ]
view(STAD_mito_VST)
nrow(STAD_mito_VST)  # Check number of mitochondrial genes)
# Save filtered matrix
write.csv(STAD_mito_VST, "mRNA_data/STAD_mito_VST.csv", row.names = TRUE)



#min-max normalization of mitochondrial genes
# Min–max normalization function
minmax_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

# Apply min–max across each gene (row)
STAD_mRNA_mito_minmax <- t(apply(STAD_mito_VST, 1, minmax_norm))
view(STAD_mRNA_mito_minmax)
#download it
write.csv(STAD_mRNA_mito_minmax, "mRNA_data/STAD_mRNA_mito_minmax.csv", row.names = TRUE)



# Apply min–max across each gene (row)
STAD_mRNA_minmax <- t(apply(STAD_vst_top, 1, minmax_norm))
view(STAD_mRNA_minmax)
#download
write.csv(STAD_mRNA_minmax, "mRNA_data/STAD_mRNA_minmax.csv", row.names = TRUE)













#read rds of miRNA file
STAD_miRNA <- readRDS("RDS files/miRNA_expression_rawCount.rds")
nrow(STAD_miRNA)# Check number of genes
ncol(STAD_miRNA) # Check number of samples)

# Remove genes with zero expression in >20% of samples
keep <- rowSums(STAD_miRNA == 0) / ncol(STAD_miRNA) <= 0.2
STAD_miRNA_filtered <- STAD_miRNA[keep, ]
nrow(STAD_miRNA_filtered)  # Check number of genes after filtering)

library(DESeq2)

# Convert to matrix if not already
STAD_miRNA_mat <- as.matrix(STAD_miRNA_filtered)

# Create a dummy sample metadata (colData) as DESeq2 requires it
sample_info <- data.frame(
  row.names = colnames(STAD_miRNA_mat),
  condition = rep("A", ncol(STAD_miRNA_mat))  # placeholder, not used for VST
)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = STAD_miRNA_mat,
  colData = sample_info,
  design = ~1
)

# Apply VST
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)

# Extract transformed matrix
STAD_miRNA_vst <- assay(vsd)
#save stad
write.csv(STAD_miRNA_vst, "miRNA_data/STAD_miRNA_vst.csv", row.names = TRUE)


# View result
View(STAD_miRNA_vst)



mitomiRNAdb <- read_csv("mitomiRNAdb.csv")
View(mitomiRNAdb)

mitomiRNAdb$tissues

# Filter for human miRNAs
mito_human <- mitomiRNAdb %>% 
  filter(grepl("Human", mth_organisms, ignore.case = TRUE))

View(mito_human)

mito_human_list <- mito_human$mirbase_id
rownames(STAD_miRNA_vst)[1:10] 

head(rownames(STAD_miRNA_vst))
head(mito_human_list)

library(miRBaseConverter)
# Convert precursor names to mature miRNAs
precursor_to_mature <- miRNA_PrecursorToMature(rownames(STAD_miRNA_vst))

# The output has a data frame with PrecursorName → MatureName mapping
head(precursor_to_mature)
# Get unique mature names for your dataset
STAD_miRNA_mature_names <- unique(precursor_to_mature$Mature1)

# Subset using mitochondrial miRNA list
STAD_miRNA_mito <- STAD_miRNA_vst[STAD_miRNA_mature_names %in% mito_human_list, ]
view(STAD_miRNA_mito)

#download
write.csv(STAD_miRNA_mito, "miRNA_data/STAD_miRNA_mito.csv", row.names = TRUE)



#min-max normalization of mitochondrial genes
# Min–max normalization function
minmax_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

# Apply min–max across each gene (row)
STAD_miRNA_mito_minmax <- t(apply(STAD_miRNA_mito, 1, minmax_norm))

# Check results
View(STAD_miRNA_mito_minmax)

STAD_miRNA_vst

# Apply min–max across each gene (row)
STAD_miRNA_minmax <- t(apply(STAD_miRNA_vst, 1, minmax_norm))

# Check results
View(STAD_miRNA_minmax)
#save file
write.csv(STAD_miRNA_minmax, "miRNA_data/STAD_miRNA_minmax.csv", row.names = TRUE)
write.csv(STAD_miRNA_mito_minmax, "miRNA_data/STAD_miRNA_mito_minmax.csv", row.names = TRUE)




