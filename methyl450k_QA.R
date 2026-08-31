
#############################################
# Methylation data preprocessing: QA Analysis
#############################################

# Packages
library(sesame)       # ‘1.30.1’
library(BiocParallel) # ‘1.46.0’
library(tidyverse)    # ‘2.0.0’
library(vroom)        # ‘1.7.1’


# 1. Load DATA
## get prefixes from files
> pfxs <- searchIDATprefixes(dir.name = "/STORAGE/csbig/multiomics-Arturo/methyl_data/") # 739

## Read raw data
> idats_raw_739 <- bplapply(
                X = pfxs,
              FUN = readIDATpair,
  BPPARAM = MulticoreParam(workers = 40) # linux only
)

# 2. Filter raw data (individuos de interés, las 3 ómicas/fenotipo AD/control = 215 expected)

# Get sample names of interest:
> methyl_samples <- methyl_biospecimen %>%
filter(individualID %in% metadata_filtered_isAD$individualID) %>%
merge(., methyl_metadata, by = "specimenID") %>%
mutate(sample = paste(Sentrix_ID, "_", Sentrix_Row_Column, sep = "")) %>%
pull(sample)

# Filter samples with methyl_samples vector:
> idats_raw_211 <- idats_raw_739[names(idats_raw_739) %in% methyl_samples]

# OBS:
De los 215 esperados se encontraron 211
-Posteriormente, será necesario volver a filtrar individuos con las 3 ómicas


# 2. QA Analysis

## Bisulfite conversion

> bis_conversion <- lapply(X = idats_raw_211, FUN = bisConversionControl)

# change list output to data frame --> For plotting:
> bis_conversion <- unlist(bis_conversion)

> bis_conversion <- data.frame(sampleID = names(bis_conversion), bis_conversion = as.numeric(bis_conversion))

# Plot bisulfite conversion:
> pdf("bis_conversion_raw_211.pdf")
> ggplot(data = bis_conversion, mapping = aes(x = bis_conversion)) +
+ geom_histogram() +
+ theme_classic()
`stat_bin()` using `bins = 30`. Pick better value `binwidth`.
> dev.off()
pdf
  2

## Detection p-value

> pvalues <- bplapply(X = idats_raw_211, FUN = function(sample){pOOBAH(sample, return.pval = TRUE)}, BPPARAM = MulticoreParam(workers = 40))

# change samples (rows) to columns:
> pvalues <- do.call(cbind, pvalues)

# Count probes with average pvalue > 0.01:
> sum(rowMeans(pvalues > 0.01) != 0, na.rm = TRUE)
[1] 96363

# Count NAs:
> sum(is.na(pvalues))
[1] 422

# Average p-value per sample:
> avg_pvalue_per_sample <- colMeans(pvalues, na.rm = TRUE)

# All avg_pvalues < 0.05?:
> all(avg_pvalue_per_sample < 0.05)
[1] TRUE

# Plot average p-value per sample:

# Format for plotting (data frame):
> avg_pvalue_per_sample_df <- data.frame(sampleID = names(avg_pvalue_per_sample), pvalue = avg_pvalue_per_sample)


# Plot:
> pdf("Detection_pvalue_per_sample_raw_211.pdf")
> ggplot(data = avg_pvalue_per_sample_df, mapping = aes(x = pvalue)) +
+ geom_density() +
+ scale_fill_viridis_d() +
+ theme_classic(base_size = 14) +
+ theme(
    legend.position = "right",
    axis.line = element_line(color = "black")
  ) +
  labs(
    x = "Detection p-value",
    y = "Density (Number of samples)",
    title = "Detection p-value Distribution"
  ) +
  xlim(c(0.000,0.055)) +
  geom_vline(xintercept = 0.05,linetype = "dashed", col = "red")
> dev.off()
pdf
  2


## Dye bias

#Plot dye bias per sample (qqplots):
> pdf("qqplots_intensityRedGrn_raw_211.pdf")
> lapply(X = idats_raw_211, FUN = sesameQC_plotRedGrnQQ)
> dev.off()

## Success probes

> pdf("plotBar_raw_211.pdf")
> sesameQC_plotBar(lapply(X = idats_raw_211, FUN = sesameQC_calcStats, "detection"))
> dev.off()
pdf
  2

## Get stats per sample:

# Calc stats:
> qc_metrics_raw <- bplapply(X = idats_raw_211, FUN = sesameQC_calcStats, BPPARAM = MulticoreParam(workers = 40))

# get stats:
> stats_raw_211 <- lapply(X = qc_metrics_raw, FUN = sesameQC_getStats)

# Get stats per sample (in data frame format):
> metric_names <- names(stats_raw_211[[1]])

> qc_stats_raw_211 <- do.call(
                    rbind,
                    lapply(
                        X = stats_raw_211,
                      FUN = function(sample_stat){
                              sample_stat[][metric_names]
                            }
                    )
                  )

> qc_stats_raw_211 <- as.data.frame(qc_stats_raw_211)

# All frac_dt >= 0.95?
> all(qc_stats_raw_211$frac_dt > 0.95)
[1] FALSE

> sum(qc_stats_raw_211$frac_dt >= 0.95)
[1] 202 # Passes filter


## qqplots summary:

plot RG_ratio vs RG_distortion:
> pdf("RG_ratio_vs_RG_distortion_raw_211.pdf")
> plot(qc_stats_raw_211$RGratio, qc_stats_raw_211$RGdistort,
     xlab="RG ratio (shift)",
     ylab="RG distortion (curvature)",
     pch=16)
abline(h=0, lty=2)
abline(v=1, lty=2)
> dev.off()

## beta values


# Get betas:

> betas_raw_211 <- do.call(
  cbind,
  lapply(
    X = idats_raw_211,
    FUN = function(prefix){
      getBetas(prefix)
    }
  )
)

# Get betas in data frame format:

> betas_raw_211_df <- betas_raw_211 %>%
  as.data.frame() %>%
  rownames_to_column(var = "CpG") %>%
  pivot_longer(
    cols = -CpG,
    names_to = "sample",
    values_to = "beta"
  )


# Plot beta values (all samples):

> pdf("b_values_raw_211.pdf")
> ggplot(betas_raw_211_df, aes(x = beta, group = sample)) +
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


# 20/08/2026

## PCA m values:

# beta values to m values:
> m_values_raw_211 <- BetaValueToMValue(
  b = betas_raw_211
)

# Filter probes with NAs, infinite, NAs:
> probes_to_keep <- <- apply(
  X = m_values_raw_211,
  MARGIN =  1,
  FUN =  function(x){
    all(is.finite(x))
  }

# Filter probes:
> m_values_raw_211_pca <- m_values_raw_211[probes_to_keep,]

# PCA (m values):
> pca_m_raw_211 <- prcomp(
  x = t(m_values_raw_211_pca),
  scale. = TRUE)

# Plot PCA:

# -> Change PCA format for plotting (prcomp -> data.frame):
> pca_m_raw_211_df <- data.frame(
  sample = rownames(pca_m_raw_211$x),
  X = pca_m_raw_211$x[,1],
  Y = pca_m_raw_211$x[,2]
)

# -> Filter metadata:
metadata_filtered_isAD_methyl_211 <- merge(x = metadata_filtered_isAD, y = methyl_biospecimen, by = "individualID") %>% merge( x = . , y = methyl_metadata, by = "specimenID") %>% mutate(sampleID = paste(Sentrix_ID, "_", Sentrix_Row_Column, sep = "")) %>% filter(sampleID %in% pca_samples_raw_211)

# Order metadata as pca_samples, so  when we  add colours, it matches correctly:
> metadata_filtered_isAD_methyl_211 <- metadata_filtered_isAD_methyl_211[match(x = pca_m_raw_211_df$sample, table = metadata_filtered_isAD_methyl_211$sampleID), ]

# Confirm that order is correct:
> all(metadata_filtered_isAD_methyl_211$sampleID == pca_m_raw_211_df$sample)
[1] TRUE

# Save metadata file:
#> write.table(x = metadata_filtered_isAD_methyl_211, file = "metadata_filtered_isAD_methyl_211.tsv", sep = "\t")

# Plot PCA (color = batch):

> pdf("pca_m_raw_batch_211.pdf")
> pca_m_raw_211_df %>%
 ggplot(mapping = aes(x = X, y = Y)
      ) +
      geom_point() +
      aes(colour = as.factor(metadata_filtered_isAD_methyl_211$batch)) +
      scale_color_discrete(name = "Batch") +
      xlab(paste("PC1 - ", pca_var_raw_211_per[1], "%", sep = "")) +
      ylab(paste("PC2 - ", pca_var_raw_211_per[2], "%", sep = "")) +
      theme_classic() +
      ggtitle("PCA") +
      stat_ellipse(geom = "polygon", aes(fill = as.factor(metadata_filtered_isAD_methyl_211$batch)), alpha = 0.2, show.legend = FALSE)
dev.off()
pdf
  2

# plot PCA (colour = isAD)
# plot PCA (colour = Sample_Plate)
# plot PCA (colour = Study-ROS_MAP)
# plot PCA (colour = Sentrix_ID)


# Fin QA
