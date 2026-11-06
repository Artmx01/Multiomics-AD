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

blas_set_num_threads(18)
omp_set_num_threads(18)


# Tune penalty param

penalties_permut_500 <- rgcca_permutation(
  blocks     = sgcca_data,                      # data
  ncomp      = 1,                               # 1 comp
  par_type   = "sparsity",                      # Parameter to tune
  par_length = 25,                              # 25 penalties
  n_perms    = 500,                              # default value (n = 20)
  n_cores    = 4,                               # Parallelize
  scale      = TRUE,                            # Data not scaled
  method     = "sgcca",                         # SGCCA method
  connection = design,                          # Design matrix (all blocks connected)
  scheme     = "centroid",                      # Centroid scheme allows negative correlation
  scale_block = "lambda1"                       # each block is divided by the square root of the highest eigenvalue of its empirical covariance matrix.
)


penalties_permut_500$best_params
#    methyl     rnaseq      mirna
# 0.08490206 0.09080668 0.13548068


pdf("penalties_permut_500_sgcca_data.pdf")
plot(penalties_permut_500)
dev.off()


# Tune ncomp

ncomps <- expand.grid(
            methyl = 1:3,
            rnaseq = 1:3,
            mirna  = 1:3
          )


ncomp_permut_50 <- rgcca_permutation(
  blocks      = sgcca_data,                     # data
  sparsity    = penalties_permut_500$best_params,
  par_type    = "ncomp",                        # Parameter to tune
  par_value   = ncomps,                              # evaluate max 5 comps
  n_perms     = 50,                            # default value (n = 20)
  n_cores     = 4,                              # Parallelize
  scale       = TRUE,                           # Data not scaled
  method      = "sgcca",                        # SGCCA method
  connection  = design,                         # Design matrix
  scheme      = "centroid",                     # Centroid scheme allows negative correlation. maximization of the sum of the absolute values of the covariances
  scale_block = "lambda1"                     # each block is divided by the square root of the highest eigenvalue of its empirical covariance matrix.
)




ncomp_permut_50$best_params
# methyl rnaseq  mirna
#    1      1      1

pdf("ncomp_permut_50_sgcca_data.pdf")
plot(ncomp_permut_50)
dev.off()















