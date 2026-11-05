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




ncomp_AD_permut_20$best_params
# methyl rnaseq  mirna
#    1      1      1

pdf("ncomp_AD_permut_20_sgcca_data.pdf")
plot(ncomp_AD_permut_20)
dev.off()




# fit SGCCA_MODELS

fit_sgcca_AD <- rgcca(
  blocks      = sgcca_AD_data,
  connection  = design,
  ncomp       = ncomp_AD_permut_20$best_params,
  scheme      = "centroid",
  scale       = TRUE,
  scale_block = "lambda1",
  method      = "sgcca",
  sparsity    = penalties_permut_500$best_params,
)

summary(fit_sgcca_AD)
Call: method='sgcca', superblock=FALSE, scale=TRUE, scale_block='lambda1',
init='svd', bias=TRUE, tol=1e-08, NA_method='na.ignore', ncomp=c(1,1,1),
response=NULL, comp_orth=TRUE
There are J = 3 blocks.
The design matrix is:
       methyl rnaseq mirna
methyl      0      1     1
rnaseq      1      0     1
mirna       1      1     0

The centroid scheme is used.
Sum_{j,k} c_jk g(cov(X_j a_j, X_k a_k) = 0.1917

The sparsity parameter used for methyl is: 0.08490206 (with 4453
variables selected)
The sparsity parameter used for rnaseq is: 0.09080668 (with 208
variables selected)
The sparsity parameter used for mirna is: 0.1354807 (with 8 variables selected)


                 

fit_sgcca_control <- rgcca(
  blocks      = sgcca_control_data, 
  connection  = design, 
  ncomp       = ncomp_AD_permut_20$best_params,
  scheme      = "centroid", 
  scale       = TRUE, 
  scale_block = "lambda1", 
  method      = "sgcca", 
  sparsity    = penalties_permut_500$best_params, 
)
                 
summary(fit_sgcca_control)
Call: method='sgcca', superblock=FALSE, scale=TRUE, scale_block='lambda1',
init='svd', bias=TRUE, tol=1e-08, NA_method='na.ignore', ncomp=c(1,1,1),
response=NULL, comp_orth=TRUE
There are J = 3 blocks.
The design matrix is:
       methyl rnaseq mirna
methyl      0      1     1
rnaseq      1      0     1
mirna       1      1     0

The centroid scheme is used.
Sum_{j,k} c_jk g(cov(X_j a_j, X_k a_k) = 0.1983

The sparsity parameter used for methyl is: 0.08490206 (with 4523
variables selected)
The sparsity parameter used for rnaseq is: 0.09080668 (with 206
variables selected)
The sparsity parameter used for mirna is: 0.1354807 (with 10
variables selected)















