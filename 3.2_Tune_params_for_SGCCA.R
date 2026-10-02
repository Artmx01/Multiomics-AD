###############################
# Tune Parameters for SGCCA
###############################


# libraries
library(RGCCA)        # ‘3.0.3’



#            -- Workflow --
#
# ----- 1.Tune Params for AD Data -----
# 1.1 Tune ncomp param (AD data)
# 1.2 Tune penalties param (AD data)
#
# ----- 2.Tune Params for Control Data -----
# 2.1 Tune ncomp param (Control data)
# 2.2 Tune penalties param (Contorl data)



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

# Penalities to evaluate:
# For SGCCA, The minimum penalty must be greater than 1/sqrt(number_column) (i.e., 0.0569 for block 3).
aux <- c(seq(from = 0.057, to = 0.09, by = 0.003), seq(from = 0.1, to  = 0.9, by = 0.1))

penalties <- NULL

for (i in aux) {
  penalties <- c(penalties, rep(i, 3))
}

penalties <- matrix(
        data = penalties, 
        ncol = 3, 
        byrow = TRUE
)



# ------------------------- 1.Tune Params for AD Data -------------------------

# 1.1 Tune ncomp param (AD data):

ncomp_AD_permut <- rgcca_permutation(
  blocks     = sgcca_AD_data,                  # data
  par_type   = "ncomp",                        # Parameter to tune
  par_value  = 10,                             # evaluate max 10 comps
  par_length = 10,                             # 1 by 1
  n_perms    = 50,                             # default value (n = 20)
  n_cores    = 2,                              # Parallelize
  scale      = TRUE,                           # Data not scaled
  method     = "sgcca",                        # SGCCA method
  connection = design,                         # Design matrix
  scheme     = "centroid",                     # Centroid scheme allows negative correlation
  scale_block = "lambda1"                      # each block is divided by the square root of the highest eigenvalue of its empirical covariance matrix.
)
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=16h 18m 10s

ncomp_AD_permut$best_params
# methyl rnaseq  mirna
#     3      3      3


# 1.2 Tune penalties param (AD data):
penalties_AD_permut <- rgcca_permutation(
  blocks     = sgcca_AD_data,                  # data
  ncomp      = ncomp_AD_permut$best_params,
  par_type   = "sparsity",                     # Parameter to tune
  par_value  = penalties,                      # Penalties to evaluate
  n_perms    = 50,                             # default value (n = 20)
  n_cores    = 2,                              # Parallelize
  scale      = TRUE,                           # Data not scaled
  method     = "sgcca",                        # SGCCA method
  connection = design,                         # Design matrix
  scheme     = "centroid",                     # Centroid scheme allows negative correlation
  scale_block = "lambda1"                      # each block is divided by the square root of the highest eigenvalue of its empirical covariance matrix.
)


# ------------------------- 2.Tune Params for Control Data -------------------------

# 2.1 Tune ncomp param (Control data):
ncomp_control_permut <- rgcca_permutation(
  blocks     = sgcca_control_data,             # data
  par_type   = "ncomp",                        # Parameter to tune
  par_value  = 10,                             # evaluate max 10 comps
  par_length = 10,                             # 1 by 1
  n_perms    = 50,                             # default value (n = 20)
  n_cores    = 2,                              # Parallelize
  scale      = TRUE,                           # Data not scaled
  method     = "sgcca",                        # SGCCA method
  connection = design,                         # Design matrix
  scheme     = "centroid",                     # Centroid scheme allows negative correlation
  scale_block = "lambda1"                      # each block is divided by the square root of the highest eigenvalue of its empirical covariance matrix.
)


# 2.2 Tune penalties param (Control data):
penalties_control_permut <- rgcca_permutation(
  blocks     = sgcca_control_data,             # data
  ncomp      = ncomp_control_permut$best_params,
  par_type   = "sparsity",                     # Parameter to tune
  par_value  = penalties,                      # Penalties to evaluate
  n_perms    = 50,                             # default value (n = 20)
  n_cores    = 2,                              # Parallelize
  scale      = TRUE,                           # Data not scaled
  method     = "sgcca",                        # SGCCA method
  connection = design,                         # Design matrix
  scheme     = "centroid",                     # Centroid scheme allows negative correlation
  scale_block = "lambda1"                      # each block is divided by the square root of the highest eigenvalue of its empirical covariance matrix.
)
