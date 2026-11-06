####################################
# Fit SGCCA for AD data & stability
####################################


# libraries
library(RGCCA)        # ‘3.0.3’
library(RhpcBLASctl)  # ‘0.23.42’








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


boot_sgcca_AD <- rgcca_bootstrap(
                  fit_sgcca_AD,
                  n_boot             = 1000,
                  n_cores            = 50,
                  balanced           = TRUE,
                  keep_all_variables = FALSE,
                  verbose            = TRUE
)


plot(
  boot_sgcca_AD, 
  block = 1, 
  n_mark = 4500, 
  display_order = TRUE
)


fit_stab_sgcca_AD <- rgcca_stability(
                      fit_sgcca_AD,
               keep = sapply(fit_sgcca_AD$a, function(x) mean(x != 0)),
            n_cores = 4, 
             n_boot = 1000,
            verbose = TRUE
)


# Sigue bootstrap con fit_stab:



  

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

  
fit_stab_sgcca_control <- rgcca_stability(
                      fit_sgcca_control,
               keep = sapply(fit_sgcca_control$a, function(x) mean(x != 0)),
            n_cores = 4, 
             n_boot = 1000,
            verbose = TRUE
)





