
######################################################
# Methylation data preprocessed: QA Analysis (post QC)
######################################################

# Packages
library(sesame)       # ‘1.30.1’
library(BiocParallel) # ‘1.46.0’
library(tidyverse)    # ‘2.0.0’
library(vroom)        # ‘1.7.1’


# Bisulfite conversion:
bis_conversion_processed_198 <- lapply(X = idats_processed_198, FUN = bisConversionControl)

# Success probes:
pdf("plotBar_processed_198.pdf")
sesameQC_plotBar(lapply(X = idats_processed_198, FUN = sesameQC_calcStats, "detection"))
dev.off()

# Get stats per sample:

## Calc stats:
qc_metrics_processed <- bplapply(X = idats_processed_198, FUN = sesameQC_calcStats, BPPARAM = MulticoreParam(workers = 40))

## Get stats:
stats_processed_198 <- lapply(X = qc_metrics_processed, FUN = sesameQC_getStats)

## Change format to data frame:
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


# plot RG ratio vs RG distortion:
pdf("RG_ratio_vs_RG_distortion_processed_198.pdf")
plot(qc_stats_processed_198$RGratio, qc_stats_processed_198$RGdistort,
     xlab="RG ratio (shift)",
     ylab="RG distortion (curvature)",
     pch=16)
abline(h=0, lty=2)
abline(v=1, lty=2)
dev.off()

pdf("RG_ratio_vs_RG_distortion_processed_198_rawScale.pdf")
plot(qc_stats_processed_198$RGratio, qc_stats_processed_198$RGdistort,
     xlab="RG ratio (shift)",
     ylab="RG distortion (curvature)",
     pch=16,
     xlim = c(0.6,1.6),
     ylim = c(1.1,1.4)
)
abline(h=0, lty=2)
abline(v=1, lty=2)
dev.off()

# qqplots RG bias corrected:
pdf("qqplots_intensityRedGrn_processed_198.pdf")
lapply(X = idats_processed_198, FUN = sesameQC_plotRedGrnQQ)
dev.off()

# Plot b values

# Change object to data frame:
betas_processed_198_filtered_pval_nosex_noCrossReactive_df <- betas_processed_198_filtered_pval_nosex_noCrossReactive %>%
  as.data.frame() %>%
  rownames_to_column(var = "CpG") %>%
  pivot_longer(
    cols = -CpG,
    names_to = "sample",
    values_to = "beta"
  )

# Plot b values (all samples):
pdf("b_values_processed_198.pdf")
ggplot(betas_processed_198_filtered_pval_nosex_noCrossReactive_df, aes(x = beta, group = sample)) +
  geom_density(alpha = 0.2, size = 0.3) +
  scale_fill_viridis_d() +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    axis.line = element_line(color = "black")
  ) +
  labs(
    x = "Beta Values",
    y = "Density (Number of probes)",
    title = "Distribution Beta Values"
  )
dev.off()
 


