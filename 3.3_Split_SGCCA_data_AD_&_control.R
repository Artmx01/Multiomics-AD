############################################
# Split sgcca data into AD and control data
############################################

# library
library(tidyverse)


#      -- Workflow --
#
# 1.1 Get AD and control sgcca subjects
# 1.2 Get sgcca AD data
# 1.3 Get sgcca control data


# ------------------------ 1.1 Get AD and control sgcca subjects -------------------------

# Get sgcca AD subjects
subjects_AD_sgcca <- methyl_metadata %>% 
                      filter(individualID %in% subjects_sgcca) %>%
                      filter(is_AD == "AD-NC_SYM") %>% 
                      pull(individualID)

length(subjects_AD_sgcca)
# [1] 102


# Get sgcca control subjects
subjects_control_sgcca <- methyl_metadata %>%
                            filter(individualID %in% subjects_sgcca) %>%
                            filter(is_AD == "control") %>% 
                            pull(individualID)

length(subjects_control_sgcca)
# [1] 71



# ------------------------ 1.2 Get sgcca AD data -------------------------


# Filter datasets for AD subjects
methyl_AD_data <- methyl_sgcca_data[subjects_AD_sgcca, ]
rnaseq_AD_data <- rnaseq_sgcca_data[subjects_AD_sgcca, ]
mirna_AD_data <- mirna_sgcca_data[subjects_AD_sgcca, ]

# rownames btw datasets must match:
identical(rownames(methyl_AD_data), rownames(rnaseq_AD_data))
# [1] TRUE

identical(rownames(methyl_AD_data), rownames(mirna_AD_data))
# [1] TRUE

# Buil sgcca_AD_data object
sgcca_AD_data <- list(
  methyl = methyl_AD_data, 
  rnaseq = rnaseq_AD_data,
  mirna  = mirna_AD_data
)

# Dims
lapply(sgcca_AD_data, dim)
# $methyl
# [1]    102 341452

# $rnaseq
# [1]   102 15045

# $mirna
# [1] 102 309

# 1.3 Get sgcca control data



# ------------------------ 1.3 Get sgcca control data -------------------------

# Filter datasets for control subjects
methyl_control_data <- methyl_sgcca_data[subjects_control_sgcca, ]
rnaseq_control_data <- rnaseq_sgcca_data[subjects_control_sgcca, ]
mirna_control_data <- mirna_sgcca_data[subjects_control_sgcca, ]

# rownames btw datasets must match:
identical(rownames(methyl_control_data), rownames(rnaseq_control_data))
# [1] TRUE

identical(rownames(methyl_control_data), rownames(mirna_control_data))
# [1] TRUE

# Buil sgcca_control_data object
sgcca_control_data <- list(
  methyl = methyl_control_data,
  rnaseq = rnaseq_control_data,
  mirna  = mirna_control_data
)

# Dims
lapply(sgcca_control_data, dim)
# $methyl
# [1]     71 341452

# $rnaseq
# [1]    71 15045

# $mirna
# [1]  71 309
