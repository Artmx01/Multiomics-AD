#########################################
# Build Matrices for SGCCA (AD & Control)
#########################################


# Packages
library(vroom)     # ‘1.7.1’
library(tidyverse) # ‘2.0.0’



#            -- Workflow --
#
# ----- 1. Load filtered data -----
# 1.1 Load methyl data & metadata
# 1.2 Load RNAseq data & metadata
# 1.3 Load miRNA data & metadata
#
# ---- 2. Build matrices for SGCCA ----
# 2.1 Split AD & Control data
# 2.2 Build Control matrix
# 2.2 Build AD matrix



# ------------------------- 1. Load filtered data -------------------------


# 1.1 Load methyl data & metadata:

# Methyl data
methyl_data <- readRDS(file = "/STORAGE/csbig/multiomics-Arturo/methyl_data/preprocessing/preprocessed_data/ROSMAP_methyl450k_mvalues_filtered_198.rds")

dim(methyl_data)
# [1] 341452    198

# Methyl metadata
methyl_metadata <- vroom(file = "/STORAGE/csbig/multiomics-Arturo/methyl_data/preprocessing/preprocessed_data/ROSMAP_methyl450k_metadata_filtered_198.tsv")


# 1.2 Load RNAseq data & metadata:

# RNAseq data
rnaseq_data <- readRDS(file = "/STORAGE/csbig/multiomics-Arturo/mRNA/preprocessing/ROSMAP_RNAseq_counts_filtered_215.rds")

dim(rnaseq_data)
# [1] 15045   215

# RNAseq metadata
rnaseq_metadata <- vroom( file = "/STORAGE/csbig/multiomics-Arturo/mRNA/preprocessing/ROSMAP_RNAseq_metadata_filtered_215.txt")


# 1.3 Load miRNA data & metadata:

# miRNA metadata
mirna_data <- readRDS("/STORAGE/csbig/multiomics-Arturo/miRNA_data/preprocessing/ROSMAP_mirna_counts_filtered_215.rds")

dim(mirna_data)
# [1] 309 188

# miRNA metadata
mirna_metadata <- vroom(file = "/STORAGE/csbig/multiomics-Arturo/miRNA_data/preprocessing/ROSMAP_mirna_metadata_filtered_215.txt")



# ------------------------- 2. Build matrices for SGCCA -------------------------


# Change methyl data colnames (sampleID) to individualID:

identical(colnames(methyl_data), methyl_metadata$sampleID)
# [1] TRUE

methyl_individualIDs <- methyl_metadata %>% pull(individualID)

colnames(methyl_data) <- methyl_individualIDs


# Change RNAseq data colnames (specimenID) to individualID:

identical(colnames(rnaseq_data), rnaseq_metadata$specimenID)
# [1] TRUE

rnaseq_individualIDs <- rnaseq_metadata %>% pull(individualID)

colnames(rnaseq_data) <- rnaseq_individualIDs


# Change miRNA data colnames (mirna_id...39) to individual ID:

identical(colnames(mirna_data), mirna_metadata$mirna_id...39)
# [1] TRUE

mirna_individualIDs <- mirna_metadata %>% pull(individualID)

colnames(mirna_data) <- mirna_individualIDs


# Get data with 3 omics:
subjects_sgcca <- intersect(colnames(methyl_data), colnames(rnaseq_data)) %>% intersect(., colnames(mirna_data))

# length(subjects_sgcca)
# [1] 173
# From which:
# 102 AD
# 71 noAD

# Filter data with 3 omics
methyl_sgcca_data <- methyl_data[,subjects_sgcca]

rnaseq_sgcca_data <- rnaseq_data[, subjects_sgcca]

mirna_sgcca_data <- mirna_data[, subjects_sgcca]

# dim(methyl_sgcca_data)
# [1] 341452    173
#
# dim(rnaseq_sgcca_data)
# [1] 15045   173
#
# dim(mirna_sgcca_data)
# [1] 309 173


# Transpose data
# Reason is sgcca object needs samples in rows and features in columns
methyl_sgcca_data <- t(methyl_sgcca_data)
rnaseq_sgcca_data <- t(rnaseq_sgcca_data)
mirna_sgcca_data <- t(mirna_sgcca_data)

# Order of rows between datasets MUST BE the same
identical(rownames(methyl_sgcca_data), rownames(rnaseq_sgcca_data))
# [1] TRUE
identical(rownames(methyl_sgcca_data), rownames(mirna_sgcca_data))
# [1] TRUE

# Build sgcca object input:
sgcca_data <- list(
  methyl = methyl_sgcca_data, 
  rnaseq = rnaseq_sgcca_data,
  mirna  = mirna_sgcca_data
)

lapply(sgcca_data, dim)
# $methyl
# [1]    173 341452

# $rnaseq
# [1]   173 15045

# $mirna
# [1] 173 309

