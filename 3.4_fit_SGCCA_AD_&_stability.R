####################################
# Fit SGCCA for AD data & stability
####################################


# libraries
library(RGCCA)        # ‘3.0.3’
library(RhpcBLASctl)  # ‘0.23.42’
library(tidyverse)
library(parallel)     # ‘4.6.0’








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

blas_set_num_threads(6)
omp_set_num_threads(6)

# Make cluster
num_cores <- 5
cl <- makeCluster(num_cores)

# Export vars into nodes
clusterExport(cl, varlist = c("sgcca_AD_data", "design", "ncomp_permut_50", "penalties_permut_500"))

# Load RGCCA package to each node
clusterEvalQ(cl, library(RGCCA))

# Function to parallelize
run_sgcca_subsample <- function(i) {
  
  # Generar el submuestreo
  random_subsample <- sample(
    x    = 1:nrow(sgcca_AD_data[[1]]), 
    size = (nrow(sgcca_AD_data[[1]]) / 2)
  )
  
  subsample_data <- list(
    methyl = sgcca_AD_data[["methyl"]][random_subsample, ],
    rnaseq = sgcca_AD_data[["rnaseq"]][random_subsample, ],
    mirna  = sgcca_AD_data[["mirna"]][random_subsample, ]
  )
  
  # Check rownames between blocks match
  if (
    !identical(rownames(subsample_data[["methyl"]]), rownames(subsample_data[["rnaseq"]])) ||
    !identical(rownames(subsample_data[["methyl"]]), rownames(subsample_data[["mirna"]]))
  ) {
    stop("rownames between blocks do not match")
  }
  
  # Run SGCCA model (with random subsample, 50%)
  sgcca_subsample <- rgcca(
    blocks      = subsample_data,
    connection  = design,
    ncomp       = ncomp_permut_50$best_params,
    sparsity    = penalties_permut_500$best_params,
    scheme      = "centroid",
    scale       = TRUE,
    scale_block = "lambda1",
    method      = "sgcca"
  )
  
  # Get list of block weights
  return(
    
    list(
    methyl = sgcca_subsample$a$methyl,
    rnaseq = sgcca_subsample$a$rnaseq,
    mirna  = sgcca_subsample$a$mirna
    )
  
  )
}

# Paralellize runs
block_weights <- parLapply(cl, 1:100, run_sgcca_subsample)

stopCluster(cl)


# Get weights per omic block
methyl_weights_AD <- do.call(cbind, lapply(block_weights, function(data) data.frame(weight = data$methyl)))
rnaseq_weights_AD <- do.call(cbind, lapply(block_weights, function(data) data.frame(weight = data$rnaseq)))
mirna_weights_AD  <- do.call(cbind, lapply(block_weights, function(data) data.frame(weight = data$mirna)))













  

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

                                           
# Make cluster
num_cores <- 5
cl <- makeCluster(num_cores)

# Export vars into nodes
clusterExport(cl, varlist = c("sgcca_control_data", "design", "ncomp_permut_50", "penalties_permut_500"))

# Load RGCCA package to each node
clusterEvalQ(cl, library(RGCCA))

# Function to parallelize
run_sgcca_subsample <- function(i) {
  
  # Generar el submuestreo
  random_subsample <- sample(
    x    = 1:nrow(sgcca_control_data[[1]]), 
    size = (nrow(sgcca_control_data[[1]]) / 2)
  )
  
  subsample_data <- list(
    methyl = sgcca_control_data[["methyl"]][random_subsample, ],
    rnaseq = sgcca_control_data[["rnaseq"]][random_subsample, ],
    mirna  = sgcca_control_data[["mirna"]][random_subsample, ]
  )
  
  # Check rownames between blocks match
  if (
    !identical(rownames(subsample_data[["methyl"]]), rownames(subsample_data[["rnaseq"]])) ||
    !identical(rownames(subsample_data[["methyl"]]), rownames(subsample_data[["mirna"]]))
  ) {
    stop("rownames between blocks do not match")
  }
  
  # Run SGCCA model (with random subsample, 50%)
  sgcca_subsample <- rgcca(
    blocks      = subsample_data,
    connection  = design,
    ncomp       = ncomp_permut_50$best_params,
    sparsity    = penalties_permut_500$best_params,
    scheme      = "centroid",
    scale       = TRUE,
    scale_block = "lambda1",
    method      = "sgcca"
  )
  
  # Get list of block weights
  return(
    
    list(
    methyl = sgcca_subsample$a$methyl,
    rnaseq = sgcca_subsample$a$rnaseq,
    mirna  = sgcca_subsample$a$mirna
    )
  
  )
}

# Paralellize runs
block_weights_control <- parLapply(cl, 1:100, run_sgcca_subsample)

stopCluster(cl)


# Get weights per omic block
methyl_weights_control <- do.call(cbind, lapply(block_weights_control, function(data) data.frame(weight = data$methyl)))
rnaseq_weights_control <- do.call(cbind, lapply(block_weights_control, function(data) data.frame(weight = data$rnaseq)))
mirna_weights_contorl  <- do.call(cbind, lapply(block_weights_control, function(data) data.frame(weight = data$mirna)))








                                           







