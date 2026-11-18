####################################
# Fit SGCCA for AD data & stability
####################################


# libraries
library(RGCCA)        # ‘3.0.3’
library(RhpcBLASctl)  # ‘0.23.42’
library(tidyverse)    # ‘2.0.0’









# Fit SGCCA model for AD dataset:
fit_sgcca_AD <- rgcca(
  blocks      = sgcca_AD_data,
  connection  = design,
  ncomp       = ncomp_permut_50$best_params,
  sparsity    = penalties_permut_500$best_params,
  scheme      = "centroid",
  scale       = TRUE,
  scale_block = "lambda1",
  method      = "sgcca"
)

print(fit_sgcca_AD)
# Fitted SGCCA model.
# The algorithm converged to a stationnary point after 41 iterations.

summary(fit_sgcca_AD)
# Call: method='sgcca', superblock=FALSE, scale=TRUE, scale_block='lambda1',
# init='svd', bias=TRUE, tol=1e-08, NA_method='na.ignore', ncomp=c(1,1,1),
# response=NULL, comp_orth=TRUE
# There are J = 3 blocks.
# The design matrix is:
#       methyl rnaseq mirna
# methyl      0      1     1
# rnaseq      1      0     1
# mirna       1      1     0

# The centroid scheme is used.
# Sum_{j,k} c_jk g(cov(X_j a_j, X_k a_k) = 0.1917

# The sparsity parameter used for methyl is: 0.08490206 (with 4453
# variables selected)
# The sparsity parameter used for rnaseq is: 0.09080668 (with 208
# variables selected)
# The sparsity parameter used for mirna is: 0.1354807 (with 8 variables selected)


# Model explains > 50% variance between blocks (AVE_inner):
fit_sgcca_AD$AVE
# $AVE_X
# $AVE_X$methyl
# [1] 0.01160635

# $AVE_X$rnaseq
# [1] 0.2045678

# $AVE_X$mirna
# [1] 0.03412276


# $AVE_X_cor
# $AVE_X_cor$methyl
# [1] 0.01160635

# $AVE_X_cor$rnaseq
# [1] 0.2045678

# $AVE_X_cor$mirna
# [1] 0.03412276


# $AVE_outer
# [1] 0.02500157

# $AVE_inner
# [1] 0.5568738



# Select most stable variables

boot_sgcca_AD <- rgcca_bootstrap(
                  fit_sgcca_AD,
                  n_boot             = 1000,
                  n_cores            = 50,
                  balanced           = TRUE,
                  keep_all_variables = FALSE,
                  verbose            = TRUE
)




features_freq_AD_df <-  boot_sgcca_AD$bootstrap %>% 
                            filter(type == "weights") %>%
                            group_by(block, var) %>%
                            summarise( 
                               freq   = mean(value != 0),
                              .groups = "drop"
                            )
                            
stable_features_AD_df <- boot_sgcca_AD$stats %>%
                           filter(type == "weights") %>%
                           left_join(features_freq_AD_df, by = c("block", "var")) %>%
                           group_by(block) %>%
                           mutate (pval.adjust = p.adjust(pval, method = "BH")) %>%
                           ungroup() %>%
                           filter(freq >= 0.8, pval.adjust < 0.05)

dim(stable_features_AD_df)
# [1] 4663   16 // 4663 features (CpGs, mRNAs, miRNAs) from which:

table(stable_features_AD_df$block)
# methyl  mirna rnaseq
#  4447      8    208











  

# EN OTRO SCRIPT:


# Fit SGCCA model for control dataset:
fit_sgcca_control <- rgcca(
  blocks      = sgcca_control_data,
  connection  = design,
  ncomp       = ncomp_permut_50$best_params,
  sparsity    = penalties_permut_500$best_params,
  scheme      = "centroid",
  scale       = TRUE,
  scale_block = "lambda1",
  method      = "sgcca"
)

print(fit_sgcca_control)
# Fitted SGCCA model.
# The algorithm converged to a stationnary point after 31 iterations.

summary(fit_sgcca_control)
# Call: method='sgcca', superblock=FALSE, scale=TRUE, scale_block='lambda1',
# init='svd', bias=TRUE, tol=1e-08, NA_method='na.ignore', ncomp=c(1,1,1),
# response=NULL, comp_orth=TRUE
# There are J = 3 blocks.
# The design matrix is:
#       methyl rnaseq mirna
# methyl      0      1     1
# rnaseq      1      0     1
# mirna       1      1     0

# The centroid scheme is used.
# Sum_{j,k} c_jk g(cov(X_j a_j, X_k a_k) = 0.1983

# The sparsity parameter used for methyl is: 0.08490206 (with 4523
# variables selected)
# The sparsity parameter used for rnaseq is: 0.09080668 (with 206
# variables selected)
# The sparsity parameter used for mirna is: 0.1354807 (with 10
# variables selected)


# Model explains > 50% variance between blocks (AVE_inner):
fit_sgcca_control$AVE
# $AVE_X
# $AVE_X$methyl
# [1] 0.01670228

# $AVE_X$rnaseq
# [1] 0.2274933

# $AVE_X$mirna
# [1] 0.1201115


# $AVE_X_cor
# $AVE_X_cor$methyl
# [1] 0.01670228

# $AVE_X_cor$rnaseq
# [1] 0.2274933

# $AVE_X_cor$mirna
# [1] 0.1201115


# $AVE_outer
# [1] 0.04366314

# $AVE_inner
# [1] 0.5341872


boot_sgcca_control <- rgcca_bootstrap(
                  fit_sgcca_control,
                  n_boot             = 1000,
                  n_cores            = 50,
                  balanced           = TRUE,
                  keep_all_variables = FALSE,
                  verbose            = TRUE
)

features_freq_control_df <-  boot_sgcca_control$bootstrap %>% 
                            filter(type == "weights") %>%
                            group_by(block, var) %>%
                            summarise( 
                               freq   = mean(value != 0),
                              .groups = "drop"
                            )

stable_features_control_df <- boot_sgcca_control$stats %>%
                           filter(type == "weights") %>%
                           left_join(features_freq_control_df, by = c("block", "var")) %>%
                           group_by(block) %>%
                           mutate (pval.adjust = p.adjust(pval, method = "BH")) %>%
                           ungroup() %>%
                           filter(freq >= 0.8, pval.adjust < 0.05)


dim(stable_features_control_df)
# [1] 4508   16 // 4508 features, from which:

table(stable_features_control_df$block)
# methyl  mirna rnaseq
#  4292     10    206



# Get stable features data (AD) for network reconstruction

CpGs_stable_AD <- stable_features_AD_df %>% filter(block == "methyl") %>% pull(var)
rnaseq_stable_AD <- stable_features_AD_df %>% filter(block == "rnaseq") %>% pull(var)
mirna_stable_AD <- stable_features_AD_df %>% filter(block == "mirna") %>% pull(var)

stable_features_AD_data <- list(
                            methyl = methyl_AD_data[, CpGs_stable_AD],
                            rnaseq = rnaseq_AD_data[, rnaseq_stable_AD],
                            mirna  = mirna_AD_data[, mirna_stable_AD]
                           )



# Fit sgcca with stable features (AD)
fit_sgcca_AD_stab_features <- rgcca(
  blocks      = stable_features_AD_data,
  connection  = design,
  ncomp       = ncomp_permut_50$best_params,
  sparsity    = 1, # Keep all variables
  scheme      = "centroid",
  scale       = TRUE,
  scale_block = "lambda1",
  method      = "sgcca"
)

# Extract data blocks from fitted model:
# These blocks are preprocessed according to the values of scale/scale_block/NA_method
stable_features_AD_data_Networks <- fit_sgcca_AD_stab_features$blocks

# transpose blocks
stable_features_AD_data_Networks <- lapply(stable_features_AD_data_Networks, t)

# Save data
saveRDS(object = stable_features_AD_data_Networks, file = "stable_features_AD_data_Networks.rds")




# Get stable features data (control) for network reconstruction

CpGs_stable_control <- stable_features_control_df %>% filter(block == "methyl") %>% pull(var)
rnaseq_stable_control <- stable_features_control_df %>% filter(block == "rnaseq") %>% pull(var)
mirna_stable_control <- stable_features_control_df %>% filter(block == "mirna") %>% pull(var)

stable_features_control_data <- list(
                                  methyl = methyl_control_data[, CpGs_stable_control],
                                  rnaseq = rnaseq_control_data[, rnaseq_stable_control],
                                  mirna  = mirna_control_data[, mirna_stable_control]
                                )


# Fit sgcca with stable features (control)
fit_sgcca_control_stab_features <- rgcca(
  blocks      = stable_features_control_data,
  connection  = design,
  ncomp       = ncomp_permut_50$best_params,
  sparsity    = 1, # Keep all variables
  scheme      = "centroid",
  scale       = TRUE,
  scale_block = "lambda1",
  method      = "sgcca"
)

# Extract data blocks from fitted model:
# These blocks are preprocessed according to the values of scale/scale_block/NA_method
stable_features_control_data_Networks <- fit_sgcca_control_stab_features$blocks

# transpose blocks
stable_features_control_data_Networks <- lapply(stable_features_control_data_Networks, t)

# Save data
saveRDS(object = stable_features_control_data_Networks, file = "stable_features_control_data_Networks.rds")



# Venn Diagram of features selected by SGCCA in AD & control:

library(ggVennDiagram) # ‘1.5.7’
library(ggplot2)       # ‘4.0.3’


# Features stable in AD
stable_features_AD_names <- c(CpGs_stable_AD, rnaseq_stable_AD, mirna_stable_AD)

# Features stable in control
stable_features_control_names <- c(CpGs_stable_control, rnaseq_stable_control, mirna_stable_control)

# Venn diagram data
Venn_data <- list(AD = stable_features_AD_names, noAD = stable_features_control_names)

# Plot and save Venn diagram
pdf("VennDiagram_AD_noAD_features_selected.pdf")
ggVennDiagram(Venn_data)
dev.off()

                                      







