##########################################################################
#                 QA Analysis illumina 450k beadchip (post QC)
##########################################################################

# 13 abril 2026
# Arturo, BM

# Cargar librerias
library(sesame)
library(sesameData)
library(BiocParallel)
library(ggplot2)
library(illuminaio)
library(tidyverse)

########################
# 2. CONTROL DE CALIDAD
########################

## Sesgos de intensidades Green vs Red

# Calcular métricas para cada muestra
qc_metrics_pro<- bplapply(
  X = idats_pro_685,
  FUN = sesame::sesameQC_calcStats,
  BPPARAM = MulticoreParam(workers = 40)
) 

# Extraer métricas / Estadísticas
stats_pro <- lapply(
  X = qc_metrics_pro,
  FUN = sesameQC_getStats
) 

# Extraer nombres de las métricas
metric_names <- names(stats_pro[[1]])

# Extraer métricas para cada muestra y convertirlo a una tabla
qc_statistics_pro <- do.call(
  rbind,
  lapply(
    X = stats_pro, 
    FUN = function(sample_stat){
      sample_stat[][metric_names]
    }
  )
)

# Convertir tabla a data frame
qc_statistics_pro <- as.data.frame(qc_statistics_pro)

# Resumen qqplots
pdf("RGratio_vs_RGdistortion_pro.pdf")
plot(qc_statistics_pro$RGratio, qc_statistics_pro$RGdistort,
     xlab="RG ratio (shift)",
     ylab="RG distortion (curvature)",
     pch=16, 
     xlim = c(0.5,2.5),
     ylim = c(1, 1.5))
abline(h=0, lty=2)
abline(v=1, lty=2)
dev.off()


## Distribución de valores beta

# Convertir valores m a beta
betas_final <- MValueToBetaValue(
                  m = m_values_pro_684_noBatch
              )

# Tranformar datos para ggplot
betas_long_postQC <- betas_final %>%
  as.data.frame() %>%
  rownames_to_column(var = "CpG") %>%
  pivot_longer(
    cols = -CpG,
    names_to = "sample",
    values_to = "beta"
  )

# Graficar distribución de valores beta
pdf("b_values_postQC_results.pdf")
ggplot(betas_long_postQC, aes(x = beta, group = sample)) +
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


## PCA 

# PCA datos corregidos
pca_pro <- prcomp(
  x = t(m_values_pro_684_noBatch), 
  scale. = TRUE)

# Volver a graficar
pdf("pca_m_pro.pdf")
plot(
  x = pca_pro$x[, 1],
  y = pca_pro$x[, 2]
)
dev.off()

# Formato para ggplot
pca_m_pro_df <- data.frame(
  sample = rownames(pca_pro$x),
  X = pca_pro$x[,1],
  Y = pca_pro$x[,2]
)

pca_pro_var <- pca_pro$sdev^2

pca_pro_var_per <- round(pca_pro_var / sum(pca_pro_var) * 100, 1)

pdf("Scree_plot_pca_m_pro.pdf")
barplot(height = pca_pro_var_per, main = "Scree plot", xlab = "Principal Component", ylab = "Percent Variation")
dev.off()


# Confirmar orden entre muestras en PCA y metadata
all(rownames(pca_pro$x) == metadata_filtered_684$sampleID) # [1] TRUE
all(pca_m_pro_df$sample == metadata_filtered_684$sampleID) # [1] TRUE

# Color batch
pdf("pca_m_pro_batch.pdf")
pca_m_pro_df %>% 
  ggplot(mapping = aes(x = X, y = Y)
  ) +
  geom_point() +
  aes(colour = as.factor(metadata_filtered_684$batch)) +
  scale_color_discrete(name = "Batch") +
  xlab(paste("PC1 - ", pca_pro_var_per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca_pro_var_per[2], "%", sep = "")) +
  theme_classic() +
  ggtitle("PCA") +
  stat_ellipse(geom = "polygon", aes(fill = as.factor(metadata_filtered_684$batch)), alpha = 0.2, show.legend = FALSE)
dev.off()

# Color ceradsc
pdf("pca_m_pro_ceradsc.pdf")
pca_m_pro_df %>% 
  ggplot(mapping = aes(x = X, y = Y)
  ) +
  geom_point() +
  aes(colour = as.factor(metadata_filtered_684$ceradsc)) +
  scale_color_discrete(name = "Cerad Stage") +
  xlab(paste("PC1 - ", pca_pro_var_per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca_pro_var_per[2], "%", sep = "")) +
  theme_classic() +
  ggtitle("PCA") +
  stat_ellipse(geom = "polygon", aes(fill = as.factor(metadata_filtered_684$ceradsc)), alpha = 0.2, show.legend = FALSE)
dev.off()

# Color sample_plate
pdf("pca_m_pro_samplePlate.pdf")
pca_m_pro_df %>% 
  ggplot(mapping = aes(x = X, y = Y)
  ) +
  geom_point() +
  aes(colour = as.factor(metadata_filtered_684$Sample_Plate)) +
  scale_color_discrete(name = "Sample Plate") +
  xlab(paste("PC1 - ", pca_pro_var_per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca_pro_var_per[2], "%", sep = "")) +
  theme_classic() +
  ggtitle("PCA") +
  stat_ellipse(geom = "polygon", aes(fill = as.factor(metadata_filtered_684$Sample_Plate)), alpha = 0.2, show.legend = FALSE)
dev.off()

# Color braaksc
pdf("pca_m_pro_braaksc.pdf")
plot(
  x = pca_pro$x[, 1],
  y = pca_pro$x[, 2],
  col = as.factor(metadata_filtered_684$braaksc)
)
dev.off()

# Color msex
pdf("pca_m_pro_msex.pdf")
plot(
  x = pca_pro$x[, 1],
  y = pca_pro$x[, 2],
  col = as.factor(metadata_filtered_684$msex)
)
dev.off()

# Color study
pdf("pca_m_pro_study.pdf")
plot(
  x = pca_pro$x[, 1],
  y = pca_pro$x[, 2],
  col = as.factor(metadata_filtered_684$Study)
)
dev.off()




# Guardar datos preprocesados (posterior integración en SGCCA)
write.table(
          x = betas_final,
       file = gzfile("betas_preprocessed_684_ROSMAP.tsv.gz"),
        sep = "\t",
      quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

write.table(
  x = m_values_pro_684_noBatch,
  file = gzfile("mvalues_preprocessed_684_ROSMAP.tsv.gz"),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)


# test
betas <- vroom::vroom("mvalues_preprocessed_684_ROSMAP.tsv.gz")

vroom_write()

