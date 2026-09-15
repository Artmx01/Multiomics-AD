#############################################
# Methylation data preprocessing: QC Analysis
#############################################

library(sesame)       # ‘1.30.1’
library(BiocParallel) # ‘1.46.0’
library(tidyverse)    # ‘2.0.0’
library(vroom)        # ‘1.7.1’
library(sva)          # ‘3.60.0’
library(limma)        # ‘3.68.4’
library(kBET)         # ‘0.99.6’


#           -- Workflow --
#
# 1.--- Filtering (Sample level) -----
# 1.1 bis_conv < 0.8
# 1.2 Average detection p-value > 0.05
# 1.3 Success probe detection < 0.95 (95%)
# 1.4 Detect & Remove outliers
# 
# 2.--- Sesame Preprocessing ------
# 2.1 Sesame preprocessing
# 2.2 Get betas in matrix format
#
# 3.--- Filtering (Probe level) -----
# 3.1 Probes that have failed in one or more samples
# 3.2 Probes that fail in at least 1% of samples (considering pvalue > 0.01)
# 3.3 Probes with sex cromosome mapping
# 3.4 Cross-reactive probes
#
# 4.--- Remove Batch effect -----
# 4.1 Get beta values to m values
# 4.2 Remove known batch effect
# 4.3 Remove unknown batch effect




# 1. ----------------------- Filtering (sample level) -------------------------

# 1.1 Filter samples with bis_conversion < 0.8

range(bis_conversion$bis_conversion)
#[1] 1.034091 1.087099 
# OBS: All samples passes filter


# 1.2 Filter samples with average detection p-value > 0.05

range(avg_pvalue_per_sample_df$pvalue)
#[1] 0.003251282 0.019761264
all(avg_pvalue_per_sample_df$pvalue < 0.05)
#[1] TRUE
# OBS: All samples passes filter


# 1.3 Filter samples with success probe detection < 0.95 (95%)

success_samples <- qc_stats_raw_211 %>% filter(frac_dt >= 0.95) %>% rownames()

# Filter
idats_raw_202 <- idats_raw_211[success_samples]

length(idats_raw_202)
#[1] 202

# OBS: 202 samples passses filter, 9 do not:
#5815381002_R05C01 0.948672....
#5822038006_R04C02 0.946480....
#5822038011_R01C01 0.938333....
#5822038011_R01C02 0.949893....
#5822038011_R03C01 0.946242....
#5822054001_R01C01 0.927540....
#5822054001_R02C01 0.942683....
#5822054001_R03C02 0.947369....
#5822054001_R05C02 0.946707....
# Bad samples are in the same 4 chips: 5815381002, 5822038006, 5822038011 & 5822054001


# 1.4 Detect and remove outliers using Mahalanabis distance alghorithm
pca_scores <- pca_m_raw_211$x[,1:10]

md <- mahalanobis(
  pca_scores,
  center = colMeans(pca_scores),
  cov = cov(pca_scores)
)

threshold <- qchisq(
  0.99,
  df = ncol(pca_scores)
)

outlier_samples <- rownames(pca_scores)[md > threshold]
outlier_samples
#[1] "5815381015_R06C02" "5822038011_R01C02" "5822071001_R01C02"
#[4] "5822071001_R02C02" "6042324057_R05C02"

# Filter outlier samples:
idats_raw_198 <- idats_raw_202[!names(idats_raw_202) %in% outlier_samples]

# length(idats_raw_198)
# [1] 198
# OBS: 198 passes filter



# 2. ------------------------ Sesame Preprocessing -------------------------


# 2.1 Sesame preprocessing
# Description: 
# Mask potential bad probes ---> (qualityMask)
# Infer color channel       ---> (inferInfiniumIChannel)
# Dye bias correction       ---> (dyeBiasNL)
# pOOBAH                    ---> (Background substraction)
# Normalization             ---> (noob)

idats_processed_198 <- bplapply(
                    X = idats_raw_198,
                  FUN = function(sample){
                          noob(
                            pOOBAH(
                              dyeBiasNL(
                                inferInfiniumIChannel(
                                  qualityMask(
                                    sample
                                  )
                                )
                              )
                            )
                          )
                        },
              BPPARAM = MulticoreParam(workers = 40)
)

# 2.2 Get betas in matrix format
betas_processed_198 <- do.call(
                  cbind,
                  lapply(
                      X = idats_processed_198,
                    FUN = function(prefix){
                            getBetas(prefix)
                          }
                   )
                 )
# Matrix:
#   rows: CpG probes
#columns: samples
 


# 3. ------------------------- Filtering (probe level) -------------------------


# 3.1 Remove any probes that have failed in one or more samples

betas_processed_198_filtered <- na.omit(betas_processed_198)

dim(betas_processed_198_filtered)
# [1] 386493    198


# 3.2 Remove probes that fail in at least 1% of samples (considering pvalue > 0.01):
bad_probes <- rownames(pvalues)[rowMeans(pvalues > 0.01) >= 0.01]

## Filter bad probes:
betas_processed_198_filtered_pval <- betas_processed_198_filtered[!rownames(betas_processed_198_filtered) %in% bad_probes, ]

dim(betas_processed_198_filtered_pval)
# [1] 356321    198


# 3.3 Filter probes with sex cromosome mapping:

# Read array methylation metadata
array_meth_metadata <- vroom(file = "/STORAGE/csbig/multiomics-Arturo/methyl_data/metadata/ROSMAP_arrayMethylation_metaData.tsv")

# Get somatic probes
somatic_probes <- array_meth_metadata %>% pull(TargetID)

# Filter
betas_processed_198_filtered_pval_nosex <- betas_processed_198_filtered_pval[rownames(betas_processed_198_filtered_pval) %in% somatic_probes, ]

dim(betas_processed_198_filtered_pval_nosex)
# [1] 341995    198


# 3.4 Filter cross-reactive probes:

# Get cross reactive probes
# From:
# Zhou et al, 2016
# https://pubmed.ncbi.nlm.nih.gov/27924034/
# DOI: 10.1093/nar/gkw967

url <- "https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/HM450/archive/202209/HM450.hg19.manifest.tsv.gz"
file <- "HM450.hg19.manifest.tsv"
download.file(url = url, destfile = file)
HM450.hg19.manifest.tsv <- vroom(file = "HM450.hg19.manifest.tsv")

# Get cross-reactive probes
bad_probes <- HM450.hg19.manifest.tsv %>% filter(MASK_general == TRUE) %>% pull(probeID)

# Filter cross-reactive probes
betas_processed_198_filtered_pval_nosex_noCrossReactive <- betas_processed_198_filtered_pval_nosex[!rownames(betas_processed_198_filtered_pval_nosex) %in% bad_probes, ]

dim(betas_processed_198_filtered_pval_nosex_noCrossReactive)
# [1] 341452    198


# 4. ------------------------- Remove batch effect -------------------------

# 4.1 Get beta values to m values
m_values_processed_198 <- BetaValueToMValue(
  b = betas_processed_198_filtered_pval_nosex_noCrossReactive
)


# Filter metadata (to 198 subjects):
metadata_filtered_isAD_methyl_processed_198 <- metadata_filtered_isAD_methyl_211 %>% filter(sampleID %in% colnames(m_values_processed_198))

# Check that all samples in m_values object matches order in metadata (needed for batch effect removing):
all(metadata_filtered_isAD_methyl_processed_198$sampleID == colnames(m_values_processed_198))
# [1] TRUE

# 4.2 Remove known batch effect

# set batch
# (That's why we need same order in metadata and  samples m_values object)
batch <- metadata_filtered_isAD_methyl_processed_198$batch

## Set protecting model (this is biology, do not touch it):
model <- model.matrix(~as.factor(is_AD), metadata_filtered_isAD_methyl_processed_198)

## Remove known batch effect: batch
m_values_processed_198_noBatch <- ComBat(
                              dat = m_values_processed_198,
                            batch = batch,
                              mod = model
                            )
# Found2batches
# Adjusting for1covariate(s) or covariate level(s)
# Standardizing Data across genes
# Fitting L/S model and finding priors
# Finding parametric adjustments
# Adjusting the Data


# 4.3 Remove unknown batch effect

# set models:
model <- model.matrix(~as.factor(is_AD), metadata_filtered_isAD_methyl_processed_198)
model0 <- model.matrix(~1, data = metadata_filtered_isAD_methyl_processed_198)

# estimate hidden variation with sva:
svobj <- sva(
  dat = m_values_processed_198_noBatch,
  mod =  model,
 mod0 = model0
)
# Number of significant surrogate variables is:  22
# Iteration (out of 5 ):1  2  3  4  5 

# Remove unknown batch effect (limma)
m_values_processed_198_noBatch_unknownSVA<- removeBatchEffect(
    x = m_values_processed_198_noBatch,
  covariates = svobj$sv, design = model
)

# QC Analysis finished




##################
PCA Combat_noBatch
##################

# PCA:
pca_m_processed_198_noBatch <- prcomp(
  x = t(m_values_processed_198_noBatch),
  scale. = TRUE)

# Change format to data frame (for plotting):

pca_m_processed_198_noBatch_df <- data.frame(
  sample = rownames(pca_m_processed_198_noBatch$x),
  X = pca_m_processed_198_noBatch$x[,1],
  Y = pca_m_processed_198_noBatch$x[,2]
)

pca_var_processed_198 <- pca_m_processed_198_noBatch$sdev^2
pca_var_processed_198_per <- round(pca_var_processed_198 / sum(pca_var_processed_198) * 100, 1)

# Check order for colouring:
all(metadata_filtered_isAD_methyl_processed_198$sampleID == pca_m_processed_198_noBatch_df$sample)
#[1] TRUE


# Plot(colour = batch):

pdf("pca_m_processed_noBatch_198_colbatch.pdf")
  pca_m_processed_198_noBatch_df %>%
 ggplot(mapping = aes(x = X, y = Y)
      ) +
      geom_point() +
      aes(colour = as.factor(metadata_filtered_isAD_methyl_processed_198$batch)) +
      scale_color_discrete(name = "Batch") +
      xlab(paste("PC1 - ", pca_var_processed_198_per[1], "%", sep = "")) +
      ylab(paste("PC2 - ", pca_var_processed_198_per[2], "%", sep = "")) +
      theme_classic() +
      ggtitle("PCA") +
      stat_ellipse(geom = "polygon", aes(fill = as.factor(metadata_filtered_isAD_methyl_processed_198$batch)), alpha = 0.2, show.legend = FALSE)
dev.off()

# plot PCA (colour = isAD)
# plot PCA (colour = Sample_Plate)
# plot PCA (colour = Study-ROS_MAP)
# plot PCA (colour = Sentrix_ID)
##################################################################



#############################
PCA Combat_noBatch_unknownSVA
#############################

# PCA:
pca_m_processed_198_noBatch_unknwonSVA <- prcomp(
  x = t(m_values_processed_198_noBatch_unknownSVA),
  scale. = TRUE)

# Change format to data frame (for plotting):
pca_m_processed_198_noBatch_unknwonSVA_df <- data.frame(
  sample = rownames(pca_m_processed_198_noBatch_unknwonSVA$x),
  X = pca_m_processed_198_noBatch_unknwonSVA$x[,1],
  Y = pca_m_processed_198_noBatch_unknwonSVA$x[,2]
)


pca_var_processed_198_noBatch_unknownSVA <- pca_m_processed_198_noBatch_unknwonSVA$sdev^2
pca_var_processed_198_noBatch_unknownSVA_per <- round(pca_var_processed_198_noBatch_unknownSVA / sum(pca_var_processed_198_noBatch_unknownSVA) * 100, 1)

all(metadata_filtered_isAD_methyl_processed_198$sampleID == pca_m_processed_198_noBatch_unknwonSVA_df$sample)
#[1] TRUE

# Plot PCA (color = batch):

pdf("pca_m_processed_noBatch_unkwnonSVA_198_colbatch.pdf")
  pca_m_processed_198_noBatch_unknwonSVA_df %>%
 ggplot(mapping = aes(x = X, y = Y)
      ) +
      geom_point() +
      aes(colour = as.factor(metadata_filtered_isAD_methyl_processed_198$batch)) +
      scale_color_discrete(name = "Batch") +
      xlab(paste("PC1 - ", pca_var_processed_198_noBatch_unknownSVA_per[1], "%", sep = "")) +
      ylab(paste("PC2 - ", pca_var_processed_198_noBatch_unknownSVA_per[2], "%", sep = "")) +
      theme_classic() +
      ggtitle("PCA") +
      stat_ellipse(geom = "polygon", aes(fill = as.factor(metadata_filtered_isAD_methyl_processed_198$batch)), alpha = 0.2, show.legend = FALSE)
dev.off()


# plot PCA (colour = isAD)
# plot PCA (colour = Sample_Plate)
# plot PCA (colour = Study-ROS_MAP)
# plot PCA (colour = Sentrix_ID)


# 24/08/2026

# Quantify global structure with silhouette score:

## RAW DATA

## Check order is correct:
all(colnames(pca_m_raw_211) == metadata_filtered_isAD_methyl_211$sampleID)
#[1] TRUE

# batch:
silhouette_score_batch_m_raw_211 <- batch_sil(pca_m_raw_211, metadata_filtered_isAD_methyl_211$batch)
silhouette_score_batch_m_raw_211
#[1] 0.1677163

# AD/Control:
silhouette_score_is_AD_m_raw_211 <- batch_sil(pca_m_raw_211, as.factor(metadata_filtered_isAD_methyl_211$is_AD))
silhouette_score_is_AD_m_raw_211
#[1] -0.004783993

# Sample plate:
silhouette_score_samplePlate_m_raw_211 <- batch_sil(pca_m_raw_211, as.factor(metadata_filtered_isAD_methyl_211$Sample_Plate))
silhouette_score_samplePlate_m_raw_211
#[1] -0.1323941

# SentrixID:
silhouette_score_sentrixID_m_raw_211 <- batch_sil(pca_m_raw_211, as.factor(metadata_filtered_isAD_methyl_211$Sentrix_ID))
silhouette_score_sentrixID_m_raw_211
#[1] -0.2310685


## KNOWN BATCH

# Check order is correct:
all(metadata_filtered_isAD_methyl_processed_198$sampleID == colnames(pca_m_processed_198_noBatch))
#[1] TRUE

# Batch:
silhouette_score_batch_m_processed_198_noBatch <- batch_sil(pca_m_processed_198_noBatch, metadata_filtered_isAD_methyl_processed_198$batch)
silhouette_score_batch_m_processed_198_noBatch
#[1] 0.03846471

# AD/Control:
silhouette_score_is_AD_m_processed_198_noBatch <- batch_sil(pca_m_processed_198_noBatch, as.factor(metadata_filtered_isAD_methyl_processed_198$is_AD))
silhouette_score_is_AD_m_processed_198_noBatch
#[1] -0.001152145

# Sample plate:
silhouette_score_samplePlate_m_processed_198_noBatch <- batch_sil(pca_m_processed_198_noBatch, as.factor(metadata_filtered_isAD_methyl_processed_198$Sample_Plate))
silhouette_score_samplePlate_m_processed_198_noBatch
#[1] -0.1638075

# SentrixID
silhouette_score_sentrixID_m_processed_198_noBatch <- batch_sil(pca_m_processed_198_noBatch, as.factor(metadata_filtered_isAD_methyl_processed_198$Sentrix_ID))
silhouette_score_sentrixID_m_processed_198_noBatch
#[1] -0.4089403


## KNOWN BATCH + UNKOWN SVA(LIMMA):

# Check order is correct:
all(metadata_filtered_isAD_methyl_processed_198$sampleID == colnames(pca_m_processed_198_noBatch_unknwonSVA))
#[1] TRUE

# batch:
silhouette_score_batch_m_processed_198_noBatch_unkownSVA <- batch_sil(pca_m_processed_198_noBatch_unknwonSVA, metadata_filtered_isAD_methyl_processed_198$batch)
silhouette_score_batch_m_processed_198_noBatch_unkownSVA
#[1] 0.002686696

# AD/Control:
silhouette_score_is_AD_m_processed_198_noBatch_unkownSVA <- batch_sil(pca_m_processed_198_noBatch_unknwonSVA, as.factor(metadata_filtered_isAD_methyl_processed_198$is_AD))
silhouette_score_is_AD_m_processed_198_noBatch_unkownSVA
#[1] 0.2126924

# Sample Plate:
silhouette_score_samplePlate_m_processed_198_noBatch_unkownSVA <- batch_sil(pca_m_processed_198_noBatch_unknwonSVA, as.factor(metadata_filtered_isAD_methyl_processed_198$Sample_Plate))
silhouette_score_samplePlate_m_processed_198_noBatch_unkownSVA
#[1] -0.09693983

# SentrixID:
silhouette_score_sentrixID_m_processed_198_noBatch_unkownSVA <- batch_sil(pca_m_processed_198_noBatch_unknwonSVA, as.factor(metadata_filtered_isAD_methyl_processed_198$Sentrix_ID))
silhouette_score_sentrixID_m_processed_198_noBatch_unkownSVA
#[1] -0.4271414



