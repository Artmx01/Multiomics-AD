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

ncomp_AD_permut_20 <- rgcca_permutation(
  blocks      = sgcca_data,                     # data
  sparsity    = penalties_permut_500$best_params,
  par_type    = "ncomp",                        # Parameter to tune
  par_value   = 3,                              # evaluate max 5 comps
  par_length  = 3,                              # 1 by 1
  n_perms     = 20,                            # default value (n = 20)
  n_cores     = 4,                              # Parallelize
  scale       = TRUE,                           # Data not scaled
  method      = "sgcca",                        # SGCCA method
  connection  = design,                         # Design matrix
  scheme      = "centroid",                     # Centroid scheme allows negative correlation. maximization of the sum of the absolute values of the covariances
  scale_block = "lambda1",                     # each block is divided by the square root of the highest eigenvalue of its empirical covariance matrix.
  tol         = 1e-3
)


ncomp_AD_permut_20$best_params
# methyl rnaseq  mirna
#    1      1      1

pdf("ncomp_AD_permut_20_sgcca_data.pdf")
plot(ncomp_AD_permut_20)
dev.off()


# En otro script:

library(tidyverse)

subjects_AD_sgcca <- methyl_metadata %>% filter(individualID %in% subjects_sgcca) %>% filter(is_AD == "AD-NC_SYM") %>% pull(individualID)

length(subjects_AD_sgcca)
# [1] 102


subjects_control_sgcca <- methyl_metadata %>% filter(individualID %in% subjects_sgcca) %>% filter(is_AD == "control") %>% pull(individualID)

length(subjects_control_sgcca)
# [1] 71


# Split AD & control Data:


methyl_AD_data <- methyl_sgcca_data[subjects_AD_sgcca, ]
rnaseq_AD_data <- rnaseq_sgcca_data[subjects_AD_sgcca, ]
mirna_AD_data <- mirna_sgcca_data[subjects_AD_sgcca, ]

# Create sgcca objects for AD & control data:

identical(rownames(methyl_AD_data), rownames(rnaseq_AD_data))
[1] TRUE

identical(rownames(methyl_AD_data), rownames(mirna_AD_data))
[1] TRUE

sgcca_AD_data <- list(
  methyl = methyl_AD_data, 
  rnaseq = rnaseq_AD_data,
  mirna  = mirna_AD_data
)

lapply(sgcca_AD_data, dim)
$methyl
[1]    102 341452

$rnaseq
[1]   102 15045

$mirna
[1] 102 309


# Control data

methyl_control_data <- methyl_sgcca_data[subjects_control_sgcca, ]
rnaseq_control_data <- rnaseq_sgcca_data[subjects_control_sgcca, ]
mirna_control_data <- mirna_sgcca_data[subjects_control_sgcca, ]


identical(rownames(methyl_control_data), rownames(rnaseq_control_data))
[1] TRUE


identical(rownames(methyl_control_data), rownames(mirna_control_data))
[1] TRUE

sgcca_control_data <- list(
  methyl = methyl_control_data,
  rnaseq = rnaseq_control_data,
  mirna  = mirna_control_data
)


lapply(sgcca_control_data, dim)
$methyl
[1]     71 341452

$rnaseq
[1]    71 15045

$mirna
[1]  71 309



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















