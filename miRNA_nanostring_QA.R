#############################################
# MiRNA data preprocessing: QA Analysis
#############################################

library(vroom)       # ‘1.7.1’
library(tidyverse)   # ‘2.0.0’
library(CePa)        # ‘0.8.2’


# 1. Load metadata:

# Load mirna_assay_nanostring_metadata:

mirna_assay_nanostring_metadata <- vroom(file = "/STORAGE/csbig/multiomics-Arturo/miRNA_data/ROSMAP_assay_miRNAarray_nanostring_metadata.csv")

# 1.1. Filter metadata mirna_assay_nanostring_metadata:

subjects <- metadata_filtered_isAD %>% pull(individualID)

subjects <- mirna_biospecimen %>% filter(individualID %in% subjects) %>% pull(specimenID)

mirna_assay_nanostring_metadata_filtered_189 <- mirna_assay_nanostring_metadata %>%
        filter(specimenID %in% subjects) %>%
        filter(!is.na(mirna_id...2))

# dim(mirna_assay_nanostring_metadata_filtered_189)
# [1] 189  14
# OBS: Not all subjects match metadata (215 expected)


# 2. Load data (gct file):
mirna_processed <- read.gct(file = "/STORAGE/csbig/multiomics-Arturo/miRNA_data/ROSMAP_arraymiRNA.gct")

# dim(mirna_processed)
# [1] 309 702
# 309 miRNAs, 702 subjects

# 2.1 Filter data (only subjects in mirna_assay_nanostring_metadata_filtered_189):

# Get subjects ID:
sampleID_DLPFC <- mirna_assay_nanostring_metadata_filtered_189 %>%
                pull(mirna_id...2)

# Filter matrix data:
mirna_processed_188 <- mirna_processed[,colnames(mirna_processed) %in% sampleID_DLPFC]

# dim(mirna_processed_188)
# [1] 309 188
# OBS: We need to filter metadata again (in PCA)


#############
# QA Analysis
#############


# Range:
range(mirna_processed_188)
# [1]  2.444961 16.501511

# PCA
pca_mirna_processed_188 <- prcomp(
             x = t(mirna_processed_188), 
        scale. = TRUE
)

# PCA to data frame (PCS 1 & 2):
pca_mirna_processed_188_df <- data.frame(
  sample = rownames(pca_mirna_processed_188$x),
  X = pca_mirna_processed_188$x[,1],
  Y = pca_mirna_processed_188$x[,2]
)


# Get % variation:
pca_var_mirna_processed_188 <- pca_mirna_processed_188$sdev^2

pca_var_mirna_processed_188_per <- round(pca_var_mirna_processed_188 / sum(pca_var_mirna_processed_188) * 100 , 1)

# head(pca_var_mirna_processed_188_per)
# [1] 14.8  7.7  5.5  4.5  3.5  3.1


# Plot PCA:

# Filter metadata:

# Build metadata_all:
mirna_metadata_all_isAD_188 <- merge(x = metadata_filtered_isAD, y = mirna_biospecimen, by = "individualID") %>%
        merge( x = ., y = mirna_assay_nanostring_metadata_filtered_189, by = "specimenID") %>% 
        filter(mirna_id...2 %in% pca_mirna_processed_188_df$sample)

# Order metadata as in PCA:
mirna_metadata_all_isAD_188 <- mirna_metadata_all_isAD_188[match(x = pca_mirna_processed_188_df$sample, table = mirna_metadata_all_isAD_188$mirna_id...2), ]

# Check order:
all(mirna_metadata_all_isAD_188$mirna_id...2 == pca_mirna_processed_188_df$sample)
#[1] TRUE

pdf("pca_mirna_processed_188_isAD.pdf")
  pca_mirna_processed_188_df %>%
 ggplot(mapping = aes(x = X, y = Y)
      ) +
      geom_point() +
      aes(colour = as.factor(mirna_metadata_all_isAD_188$is_AD)) +
      scale_color_discrete(name = "AD/Control") +
      xlab(paste("PC1 - ", pca_var_mirna_processed_188_per[1], "%", sep = "")) +
      ylab(paste("PC2 - ", pca_var_mirna_processed_188_per[2], "%", sep = "")) +
      theme_classic() +
      ggtitle("PCA") +
      stat_ellipse(geom = "polygon", aes(fill = as.factor(mirna_metadata_all_isAD_188$is_AD)), alpha = 0.2, show.legend = FALSE)
dev.off()

# plot (colour = plate)
# plot (colour = msex)


