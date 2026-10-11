###############################
# Tune Parameters for SGCCA
###############################


# libraries
library(RGCCA)        # ‘3.0.3’
library(RhpcBLASctl)  # ‘0.23.42’


#            -- Workflow --
#
# ----- 1.Tune Params for AD Data -----
# 1.1 Tune ncomp param (AD data)
# 1.2 Tune penalties param (AD data)
#
# ----- 2.Tune Params for Control Data -----
# 2.1 Tune ncomp param (Control data)
# 2.2 Tune penalties param (Control data)



# Design matrix:
# With this design, All blocks are connected
design <- matrix(
  data = c(0,1,1,
           1,0,1,
           1,1,0), 
  ncol  = 3,
  nrow  = 3, 
  byrow = TRUE
)

blas_set_num_threads(10)
omp_set_num_threads(10)


# ------------------------- 1.Tune Params for AD Data -------------------------

# 1.1 Tune ncomp param (AD data):

# Previous EDA analysis (n_perms = 50; ncomp = 1,2...,10) suggested 3 components as the best param.
# Considering max 5 comps & 500 permutations

ncomp_AD_permut_500 <- rgcca_permutation(
  blocks      = sgcca_AD_data,                  # data
  par_type    = "ncomp",                        # Parameter to tune
  par_value   = 5,                              # evaluate max 5 comps
  par_length  = 5,                              # 1 by 1
  n_perms     = 500,                            # default value (n = 20)
  n_cores     = 3,                              # Parallelize
  scale       = TRUE,                           # Data not scaled
  method      = "sgcca",                        # SGCCA method
  connection  = design,                         # Design matrix
  scheme      = "centroid",                     # Centroid scheme allows negative correlation. maximization of the sum of the absolute values of the covariances
  scale_block = "lambda1"                      # each block is divided by the square root of the highest eigenvalue of its empirical covariance matrix.
)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=01d 07h 40m 59s

# ncomp_AD_permut_500$best_params
# methyl rnaseq  mirna
#     3      3      3


# 1.2 Tune penalties param (AD data):
penalties_AD_permut_500 <- rgcca_permutation(
  blocks     = sgcca_AD_data,                   # data
  ncomp      = ncomp_AD_permut_500$best_params, # Best ncomp params
  par_type   = "sparsity",                      # Parameter to tune
  par_value  = 0.3,                             # Max penalty to evaluate
  par_length = 30,                              # 30 penalties
  n_perms    = 500,                             # default value (n = 20)
  n_cores    = 3,                               # Parallelize
  scale      = TRUE,                            # Data not scaled
  method     = "sgcca",                         # SGCCA method
  connection = design,                          # Design matrix
  scheme     = "centroid",                      # Centroid scheme allows negative correlation
  scale_block = "lambda1"                       # each block is divided by the square root of the highest eigenvalue of its empirical covariance matrix.
)
saveRDS(object = penalties_AD_permut_500, file = "penalties_AD_permut_500.rds")



# ------------------------- 2.Tune Params for Control Data -------------------------

# 2.1 Tune ncomp param (Control data):

blas_set_num_threads(14)
omp_set_num_threads(14)

ncomp_control_permut_500 <- rgcca_permutation(
  blocks      = sgcca_control_data,             # data
  par_type    = "ncomp",                        # Parameter to tune
  par_value   = 5,                              # evaluate max 5 comps
  par_length  = 5,                              # 1 by 1
  n_perms     = 500,                            # default value (n = 20)
  n_cores     = 2,                              # Parallelize
  scale       = TRUE,                           # Data not scaled
  method      = "sgcca",                        # SGCCA method
  connection  = design,                         # Design matrix
  scheme      = "centroid",                     # Centroid scheme allows negative correlation
  scale_block = "lambda1"                      # each block is divided by the square root of the highest eigenvalue of its empirical covariance matrix.
)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=01d 11h 26m 34s

# ncomp_control_permut_500$best_params
# methyl rnaseq  mirna
#     1      1      1



# 2.2 Tune penalties param (Control data):
penalties_control_permut_500 <- rgcca_permutation(
  blocks     = sgcca_control_data,                   # data
  ncomp      = ncomp_control_permut_500$best_params,
  par_type   = "sparsity",                           # Parameter to tune
  par_value  = 0.3,                                  # Penalties to evaluate
  par_length = 30,
  n_perms    = 500,                                  # default value (n = 20)
  n_cores    = 2,                                    # Parallelize
  scale      = TRUE,                                 # Data not scaled
  method     = "sgcca",                              # SGCCA method
  connection = design,                               # Design matrix
  scheme     = "centroid",                           # Centroid scheme allows negative correlation
  scale_block = "lambda1"                            # each block is divided by the square root of the highest eigenvalue of its empirical covariance matrix.
)
saveRDS(object = penalties_control_permut_500, file = "penalties_control_permut_500.rds")

