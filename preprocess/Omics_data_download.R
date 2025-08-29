#miRNA expression,DNA methylation,Mutation (MAF),Copy number variation (CNV)
setwd("D:/Multiomics_STAD/STAD_multiomics")
library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)
# get a list of projects
gdcprojects <- getGDCprojects()
getProjectSummary('TCGA-STAD')



#--------------------------------
#mRNA expression
#---------------------------------
query_mRNA <- GDCquery(
  project = "TCGA-STAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor"
)
# download data - GDCdownload
GDCdownload(query_mRNA, files.per.chunk = 10, method = "api")
mRNA_data <- GDCprepare(query_mRNA, summarizedExperiment = TRUE)
STAD_matrix_normal <- assay(tcga_STAD_data_normal, 'FPKM')



# ---------------------------
# 1. miRNA expression
# ---------------------------
query_miRNA <- GDCquery(
  project = "TCGA-STAD",
  data.category = "Transcriptome Profiling",
  data.type = "miRNA Expression Quantification",
  workflow.type = "BCGSC miRNA Profiling",
  sample.type = "Primary Tumor",
)
getResults(query_miRNA)
GDCdownload(query_miRNA)
miRNA_data <- GDCprepare(query_miRNA)
str(miRNA_data)


library(reshape2)

expr_matrix <- dcast(
  miRNA_data,
  miRNA_ID ~ sample,
  value.var = "read_count"
)
view(miRNA_data)
colnames(miRNA_data)
head(miRNA_data)


# Keep only columns with RPM data
rpm_cols <- grep("^reads_per_million_miRNA_mapped", colnames(miRNA_data), value = TRUE)

# Extract matrix
expr_matrix_miRNA <- as.matrix(miRNA_data[, rpm_cols])
view(expr_matrix_miRNA)
# Set rownames to miRNA IDs
rownames(expr_matrix_miRNA) <- miRNA_data$miRNA_ID
# Optional: Clean column names to keep only sample IDs
colnames(expr_matrix_miRNA) <- sub("^reads_per_million_miRNA_mapped_", "", rpm_cols)



#save as .rds
saveRDS(expr_matrix_miRNA, file = "miRNA_expression_RPM.rds")



# ---------------------------
# 2. DNA Methylation
# ---------------------------
BiocManager::install("DMRCrate")
query_meth <- GDCquery(
  project = "TCGA-STAD",
  data.category = "DNA Methylation",
  data.type = "Methylation Beta Value",
  platform = "Illumina Human Methylation 450",
  sample.type = "Primary Tumor"
)
getResults(query_meth)

GDCdownload(query_meth, files.per.chunk = 15, method = "api")
rowData(meth_data)[1:5, ] 
# Load already downloaded data without downloading again
meth_data <- GDCprepare(query_meth, save = TRUE, save.filename = "meth_data.rda")
meth_matrix <- assay(meth_data)
rowData(meth_data)[1:5, ]
head(meth_matrix)
saveRDS(meth_matrix, file = "methylation_beta_values.rds")






## ---------------------------------------
## 1. Mutation (MAF)
## ---------------------------------------
query_mut <- GDCquery(
  project = "TCGA-STAD",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)
GDCdownload(query_mut)
maf_data <- GDCprepare(query_mut)
view(maf_data)

#save
saveRDS(maf_data, file = "mutation_maf.rds")

# Or get as MAF object for analysis
library(maftools)
maf_object <- GDCquery_Maf("STAD", pipelines = "mutect2")



## -------------------
## 2. Copy Number Variation (CNV)
## -------------------
query_cnv <- GDCquery(
  project = "TCGA-STAD",
  data.category = "Copy Number Variation",
  data.type = "Copy Number Segment"
)
GDCdownload(query_cnv)
cnv_data <- GDCprepare(query_cnv)
view(cnv_data)
saveRDS(cnv_data, file = "cnv_data.rds")





## -------------------
## 3. Proteomics (RPPA from TCPA)
## -------------------
library(TCGAbiolinks)

query_rppa <- GDCquery(
  project = "TCGA-STAD",
  data.category = "Proteome Profiling",
  data.type = "Protein Expression Quantification",  # Correct name
  sample.type = "Primary Tumor"
)

# Download
GDCdownload(query_rppa, method = "api")

# Prepare for analysis
GDCdownload(query_rppa, method = "api", files.per.chunk = 10)
rppa_data <- GDCprepare(query_rppa)
#assay
view(rppa_data)

# Check matrix
head(assay(rppa_data))

# Save as RDS
saveRDS(rppa_data, file = "rppa_data.rds")



library(TCGAbiolinks)
library(tidyverse)
# Query clinical data for STAD
clinical <- GDCquery_clinic("TCGA-STAD", type = "clinical")
head(clinical)
view(clinical)
colnames(clinical)
clinical$days_to_recurrence
#download clinical data
write.csv(clinical, "clinical_data/STAD_clinical_data.csv", row.names = FALSE)


#read clinical data
clinical <- read.csv("clinical_data/STAD_clinical_data.csv")

# read rds file
multiomics_matrix <- readRDS("multiomics_matrix/multiomics_matrix_for_AE.rds")
rownames(multiomics_matrix)[1:5]
clinical$submitter_id[1:5]


#trim multiomics matrix to 3 fields
# Example: "TCGA-BR-6457-01A" → "TCGA-BR-6457"
multiomics_ids <- substr(rownames(multiomics_matrix), 1, 12)


# Keep only clinical rows with matching IDs
clinical_filtered <- mmc1[mmc1$bcr_patient_barcode %in% multiomics_ids, ]
view(clinical_filtered)



library(readxl)
mmc1 <- read_excel("mmc1.xlsx")
View(mmc1)
colnames(clinical_filtered)


#CHECK HOW MANY PEOPLE HAVE DFI NA
sum(is.na(clinical_filtered$PFI))

#SAVE THE DATA
write.csv(clinical_filtered, "clinical_data/STAD_clinical_filtered.csv", row.names = FALSE)









