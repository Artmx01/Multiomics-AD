#########################################
# Build Matrices for SGCCA (AD & Control)
#########################################


# Packages
library(vroom)     # ‘1.7.1’
library(tidyverse) # ‘2.0.0’



#           -- Workflow --
#
# ----- 1. Load filtered data -----
# 1.1 Load methyl data & metadata
# 1.2 Load RNAseq data & metadata
# 1.3 Load miRNA data & metadata
#
# ---- 2. Build matrices for SGCCA ----
# 2.1 Split AD & Control data
# 2.1 Build Control matrix
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


# 2.1 Split AD & Control data:

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


# Split data:

# Mehtyl Control data

methyl_control_subjects <- methyl_metadata %>%
                            filter(is_AD == "control") %>%
                            pull(individualID)

methyl_control_data <- methyl_data[, methyl_control_subjects]

dim(methyl_control_data)
# [1] 341452     81


# RNAseq Control data

rnaseq_control_subjects <- rnaseq_metadata %>%
                            filter(is_AD == "control") %>%
                            pull(individualID)

rnaseq_control_data <- rnaseq_data[, rnaseq_control_subjects]

dim(rnaseq_control_data)
# [1] 15045    88


# miRNA Control data

mirna_control_subjects <- mirna_metadata %>%
                            filter(is_AD == "control") %>%
                            pull(individualID)

mirna_control_data <- mirna_data[, mirna_control_subjects]

dim(mirna_control_data)
# [1] 309  76


# Methyl AD data

methyl_AD_subjects <- methyl_metadata %>%
                        filter(is_AD == "AD-NC_SYM") %>%
                        pull(individualID)

methyl_AD_data <- methyl_data[,methyl_AD_subjects]

dim(methyl_AD_data)
# [1] 341452    117


# RNAseq AD data

rnaseq_AD_subjects <- rnaseq_metadata %>%
                        filter(is_AD == "AD-NC_SYM") %>%
                        pull(individualID)

rnaseq_AD_data <- rnaseq_data[, rnaseq_AD_subjects]

dim(rnaseq_AD_data)
# [1] 15045   127


# miRNA AD data

mirna_AD_subjects <- mirna_metadata %>%
                      filter(is_AD == "AD-NC_SYM") %>%
                      pull(individualID)

mirna_AD_data <- mirna_data[, mirna_AD_subjects]

dim(mirna_AD_data)
# [1] 309 112



# 2.1 Build Control matrix:

# Get control subjects with the 3 omics
control_subjects_sgcca <- intersect(x = colnames(methyl_control_data), y = colnames(rnaseq_control_data)) %>% 
                            intersect(x = ., y = colnames(mirna_control_data))

length(control_subjects_sgcca)
# [1] 71 // 71 Control subjects with 3 omics after data preprocessed


# Filter control datasets

methyl_control_data <- methyl_control_data[, control_subjects_sgcca]

rnaseq_control_data <- rnaseq_control_data[, control_subjects_sgcca]

mirna_control_data <- mirna_control_data[, control_subjects_sgcca]

# Transpose datasets

methyl_control_data <- t(methyl_control_data)
# [1]     71 341452

rnaseq_control_data <- t(rnaseq_control_data)
# [1]    71 15045

mirna_control_data <- t(mirna_control_data)
# [1]  71 309

# Build SGCCA control data:

# Order of rows between datasets MUST BE the same

identical(rownames(methyl_control_data), rownames(rnaseq_control_data))
#[1] TRUE

identical(rownames(methyl_control_data), rownames(mirna_control_data))
# [1] TRUE

# SGCCA control data
sgcca_control_data <- list(
  methyl = methyl_control_data, 
  rnaseq = rnaseq_control_data,
  mirna  = mirna_control_data
)


# 2.2 Build AD matrix:

















