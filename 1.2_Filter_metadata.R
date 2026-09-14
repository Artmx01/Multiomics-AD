#############################################################
# Filter metadata: Subjects with 3 omics & AD/Control (DMS-5)
#############################################################


# Packages
library(vroom)     # ‘1.7.1’
library(tidyverse) # ‘2.0.0’


# Workflow:
#----------
# 1. Identify subjects with 3 omics
# 1.1 Load metadatas
# 1.2 Load biospecimen metadatas
# 1.3 Get individual IDs per omic
# 1.4 Intersection (subjects with 3 omics)
# 2. Filter metadata with subjects selected & is_AD (DMS-5)
# 2.1 Load metadata clinic
# 2.2 Filter metadata_clinic keeping subjects with 3 omics and is_AD definition (DMS-5)
# 2.3 Save metadata filtered



# 1. Identify subjects with 3 omics:

# 1.1 Load metadatas

# methyl metadata
 methyl_metadata <- vroom(file = "/STORAGE/csbig/multiomics-Arturo/methyl_data/metadata/ROSMAP_assay_methylationArray_metadata.csv")

# Mirna metadata
 mirna_metadata <- vroom(file = "/STORAGE/csbig/multiomics-Arturo/miRNA_data/ROSMAP_assay_miRNAarray_nanostring_metadata.csv")

# rnaseq metadata
 rna_seq_metadata <- vroom(file = "/STORAGE/csbig/multiomics-Arturo/mRNA/RNA_seq_metadata_filteredQC_DLPFC.txt")


# 1.2 Load biospecimens metadata

# biospecimen metadata
metadata_biospecimen <- vroom::vroom("/STORAGE/csbig/sc_ADers/metadata/ROSMAP_biospecimen_metadata.csv")

# methyl biospecimen
methyl_biospecimen <- metadata_biospecimen %>%
   filter(assay == "methylationArray")

# mirna_biospecimen
mirna_biospecimen <- metadata_biospecimen %>%
  filter(assay == "mirnaArray")

# RNAseq biospecimen
rnaseq_biospecimen <- metadata_biospecimen %>%
 filter(assay == "rnaSeq")


# 1.3 Get individual IDs per omic

# methylation
methyl_individualID <- merge(x = methyl_metadata, y = methyl_biospecimen, by = "specimenID") %>%
  pull(individualID)

# length(methyl_individualID)
# [1] 748

# mirna
mirna_individualID <- merge(x = mirna_metadata, y = mirna_biospecimen, by = "specimenID") %>% pull(individualID)

# length(mirna_individualID)
# [1] 748

# rnaseq
rnaseq_individualID <- rna_seq_metadata %>% pull(individualID)

# length(rnaseq_individualID)
# [1] 774


# 1.4 Intersection (subjects with 3 omics)

# methyl-mirna
methyl_mirna_ID <- intersect(methyl_individualID, mirna_individualID)

#methyl-mirna-rnaseq
methyl_mirna_rnaseq_ID <- intersect(methyl_mirna_ID, rnaseq_individualID)

# length(methyl_mirna_rnaseq_ID)
# [1] 414
# OBS: 414 subjects has the 3 data omics (methylation 450k array, miRNA nanostring ncounter, RNAseq)



# 2. Filter metadata with subjects selected & is_AD (DMS-5):

# 2.1 Load metadata clinic
metadata_clinic <- vroom("/STORAGE/csbig/sc_ADers/metadata/ROSMAP_clinical.csv")

# 2.2 Filter metadata_clinic keeping subjects with 3 omics and is_AD definition (DMS-5):
metadata_filtered_isAD <- metadata_clinic %>%
  filter(individualID %in% methyl_mirna_rnaseq_ID) %>%
  mutate(is_AD = case_when(
    cogdx == 1 & ceradsc %in% c(3, 4) ~ "control",
    cogdx %in% c(4, 5) & braaksc >= 3 & ceradsc %in% c(1, 2) ~ "AD-NC_SYM",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(is_AD))

# dim(metadata_filtered_isAD)
# [1] 215  19

# table(metadata_filtered_isAD$is_AD)
# AD-NC_SYM   control
#      127        88

# Summary:
# For subsequent analysis, 215 subjects are kept, from which:
# 127 AD
# 88 control

# 2.3 Save metadata filtered
saveRDS(object = metadata_filtered_isAD, file = "metadata_filtered_multiomics_isAD_215.rds")


