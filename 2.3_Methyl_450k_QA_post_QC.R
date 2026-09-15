######################################################
# Methylation data preprocessed: QA Analysis (post QC)
######################################################


# Packages
library(sesame)       # ‘1.30.1’
library(BiocParallel) # ‘1.46.0’
library(tidyverse)    # ‘2.0.0’
library(vroom)        # ‘1.7.1’



#         -- Workflow --
#
# 1. ----- QA Analysis -----
# 1.1 Success probes
# 1.2 Get stats per sample
# 1.3 Dye bias
# 1.4 b values distribution
# 1.5 M values distribution
# 1.6 PCA



# ------------------------- QA Analysis -------------------------

# 1.1 Success probes:

pdf("plotBar_processed_198.pdf")
sesameQC_plotBar(
  lapply(
    X   = idats_processed_198, 
    FUN = sesameQC_calcStats, "detection")
)
dev.off()


# 1.2 Get stats per sample:

# Calc stats
qc_metrics_processed <- bplapply(
  X       = idats_processed_198, 
  FUN     = sesameQC_calcStats, 
  BPPARAM = MulticoreParam(workers = 40)
)

# Get stats
stats_processed_198 <- lapply(
  X   = qc_metrics_processed, 
  FUN = sesameQC_getStats
)

# Change format to data frame
metric_names <- names(stats_processed_198[[1]])

qc_stats_processed_198 <- do.call(
                    rbind,
                    lapply(
                        X = stats_processed_198,
                      FUN = function(sample_stat){
                              sample_stat[][metric_names]
                            }
                    )
                  )

qc_stats_processed_198 <- as.data.frame(qc_stats_processed_198)


# 1.3 Dye bias:

# plot RG ratio vs RG distortion (Summary of qqplots)
pdf("RG_ratio_vs_RG_distortion_processed_198.pdf")
plot(qc_stats_processed_198$RGratio, qc_stats_processed_198$RGdistort,
     xlab = "RG ratio (shift)",
     ylab = "RG distortion (curvature)",
     pch  = 16)
abline(h = 0, lty = 2)
abline(v = 1, lty = 2)
dev.off()

# plot RG ratio vs RG distortion (raw scale)
pdf("RG_ratio_vs_RG_distortion_processed_198_rawScale.pdf")
plot(qc_stats_processed_198$RGratio, qc_stats_processed_198$RGdistort,
     xlab = "RG ratio (shift)",
     ylab = "RG distortion (curvature)",
     pch  = 16,
     xlim = c(0.6,1.6),
     ylim = c(1.1,1.4)
)
abline(h = 0, lty = 2)
abline(v = 1, lty = 2)
dev.off()

# qqplots RG bias corrected
pdf("qqplots_intensityRedGrn_processed_198.pdf")
lapply(
  X   = idats_processed_198, 
  FUN = sesameQC_plotRedGrnQQ
)
dev.off()


# 1.4 b values distribution:

# Change object to data frame
betas_processed_198_filtered_pval_nosex_noCrossReactive_df <- betas_processed_198_filtered_pval_nosex_noCrossReactive %>%
  as.data.frame() %>%
  rownames_to_column(var = "CpG") %>%
  pivot_longer(
    cols      = -CpG,
    names_to  = "sample",
    values_to = "beta"
  )

# Plot b values (all samples)
pdf("b_values_processed_198.pdf")
ggplot(
  betas_processed_198_filtered_pval_nosex_noCrossReactive_df, 
  aes(x = beta, group = sample)
) +
  geom_density(alpha = 0.2, size = 0.3) +
  scale_fill_viridis_d() +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    axis.line = element_line(color = "black")
  ) +
  labs(
    x     = "Beta Values",
    y     = "Density (Number of probes)",
    title = "Distribution Beta Values"
  )
dev.off()


# 1.5 M values distribution:

# Plot m-values
pdf("m_values_processed_198.pdf")
plot(
  x    = density(m_values_processed_198_noBatch_unknownSVA),
  main = "M values distribution", 
  xlab = "m values", 
  ylab = "Density"
)
dev.off()

# Range of m values
range(m_values_processed_198_noBatch_unknownSVA)
# [1] -9.452739  8.130852


# 1.6 PCA:

#  get PCA
pca_m_processed_198_noBatch_unknwonSVA <- prcomp(
  x      = t(m_values_processed_198_noBatch_unknownSVA),
  scale. = TRUE
)

# Change format to data frame (for plotting)
pca_m_processed_198_noBatch_unknwonSVA_df <- data.frame(
  sample = rownames(pca_m_processed_198_noBatch_unknwonSVA$x),
  X = pca_m_processed_198_noBatch_unknwonSVA$x[,1],
  Y = pca_m_processed_198_noBatch_unknwonSVA$x[,2]
)

# Get % variance per PC
pca_var_processed_198_noBatch_unknownSVA <- pca_m_processed_198_noBatch_unknwonSVA$sdev^2

pca_var_processed_198_noBatch_unknownSVA_per <- round(pca_var_processed_198_noBatch_unknownSVA / sum(pca_var_processed_198_noBatch_unknownSVA) * 100, 1)

# Order is important when colouring
all(metadata_filtered_isAD_methyl_processed_198$sampleID == pca_m_processed_198_noBatch_unknwonSVA_df$sample)
# [1] TRUE

# Plot PCA (color = batch)
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
      stat_ellipse(
        geom        = "polygon", 
        aes(fill = as.factor(metadata_filtered_isAD_methyl_processed_198$batch)),
        alpha       = 0.2, 
        show.legend = FALSE
      )
dev.off()

# Plot PCA (color = is_AD)
pdf("pca_m_processed_noBatch_unkwnonSVA_198_isAD.pdf")
  pca_m_processed_198_noBatch_unknwonSVA_df %>%
  ggplot(mapping = aes(x = X, y = Y)
      ) +
      geom_point() +
      aes(colour = as.factor(metadata_filtered_isAD_methyl_processed_198$is_AD)) +
      scale_color_discrete(name = "AD/Control") +
      xlab(paste("PC1 - ", pca_var_processed_198_noBatch_unknownSVA_per[1], "%", sep = "")) +
      ylab(paste("PC2 - ", pca_var_processed_198_noBatch_unknownSVA_per[2], "%", sep = "")) +
      theme_classic() +
      ggtitle("PCA") +
      stat_ellipse(
        geom        = "polygon",
        aes(fill = as.factor(metadata_filtered_isAD_methyl_processed_198$is_AD)),
        alpha       = 0.2,
        show.legend = FALSE
      )
dev.off()

# plot PCA (colour = Sample_Plate)
pdf("pca_m_processed_noBatch_unkwnonSVA_198_SamplePlate.pdf")
  pca_m_processed_198_noBatch_unknwonSVA_df %>%
  ggplot(mapping = aes(x = X, y = Y)
      ) +
      geom_point() +
      aes(colour = as.factor(metadata_filtered_isAD_methyl_processed_198$Sample_Plate)) +
      scale_color_discrete(name = "Sample Plate") +
      xlab(paste("PC1 - ", pca_var_processed_198_noBatch_unknownSVA_per[1], "%", sep = "")) +
      ylab(paste("PC2 - ", pca_var_processed_198_noBatch_unknownSVA_per[2], "%", sep = "")) +
      theme_classic() +
      ggtitle("PCA") +
      stat_ellipse(
        geom        = "polygon",
        aes(fill = as.factor(metadata_filtered_isAD_methyl_processed_198$Sample_Plate)),
        alpha       = 0.2,
        show.legend = FALSE
      )
dev.off()

# plot PCA (colour = Study-ROS_MAP)
pdf("pca_m_processed_noBatch_unkwnonSVA_198_ROS_MAP.pdf")
  pca_m_processed_198_noBatch_unknwonSVA_df %>%
  ggplot(mapping = aes(x = X, y = Y)
      ) +
      geom_point() +
      aes(colour = as.factor(metadata_filtered_isAD_methyl_processed_198$Study)) +
      scale_color_discrete(name = "ROS-MAP") +
      xlab(paste("PC1 - ", pca_var_processed_198_noBatch_unknownSVA_per[1], "%", sep = "")) +
      ylab(paste("PC2 - ", pca_var_processed_198_noBatch_unknownSVA_per[2], "%", sep = "")) +
      theme_classic() +
      ggtitle("PCA") +
      stat_ellipse(
        geom        = "polygon",
        aes(fill = as.factor(metadata_filtered_isAD_methyl_processed_198$Study)),
        alpha       = 0.2,
        show.legend = FALSE
      )
dev.off()

# plot PCA (colour = Sentrix_ID)
pdf("pca_m_processed_noBatch_unkwnonSVA_198_SentrixID.pdf")
  pca_m_processed_198_noBatch_unknwonSVA_df %>%
  ggplot(mapping = aes(x = X, y = Y)
      ) +
      geom_point() +
      aes(colour = as.factor(metadata_filtered_isAD_methyl_processed_198$Sentrix_ID)) +
      scale_color_discrete(name = "Sentrix_ID") +
      theme_classic() +
      ggtitle("PCA") +
      stat_ellipse(
        geom        = "polygon",
        aes(fill = as.factor(metadata_filtered_isAD_methyl_processed_198$Sentrix_ID)),
        alpha       = 0.2,
        show.legend = FALSE
      )
dev.off()

# ------------------------- QA Analysis finished (post QC) -------------------------

# Save data
saveRDS(object = m_values_processed_198_noBatch_unknownSVA, file = "ROSMAP_methyl450k_mvalues_filtered_198.rds")

# Save metadata
vroom_write(x = metadata_filtered_isAD_methyl_processed_198, file = "ROSMAP_methyl450k_metadata_filtered_198.tsv")
