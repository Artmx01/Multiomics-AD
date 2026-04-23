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

# PCA de valores corregidos
pca_corrected_data <- prcomp(
  x = t(data_corrected),
  scale. = TRUE
)

# Transformar datos para ggplot
pca_corrected_data_df <- data.frame(
  sample = rownames(pca_corrected_data$x),
  X = pca_corrected_data$x[,1],
  Y = pca_corrected_data$x[,2]
)

pca_pro_var <- pca_corrected_data$sdev^2

pca_pro_var_per <- round(pca_pro_var / sum(pca_pro_var) * 100, 1)

all(pca_corrected_data_df$sample == metadata_all_filtered_357$sampleID) # [1] TRUE

# Color batch
pdf("pca_data_corrected_combat_sva_limma.pdf")
pca_corrected_data_df %>%
  ggplot(mapping = aes(x = X, y = Y)
  ) +
  geom_point() +
  aes(colour = as.factor(metadata_all_filtered_357$is_AD)) +
  scale_color_discrete(name = "AD") +
  xlab(paste("PC1 - ", pca_pro_var_per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca_pro_var_per[2], "%", sep = "")) +
  theme_classic() +
  ggtitle("PCA") +
  stat_ellipse(geom = "polygon", aes(fill = as.factor(metadata_all_filtered_357$is_AD)), alpha = 0.2, show.legend = FALSE)
dev.off()


##################################
# Evaluar si hubo sobrecorrección / Cuantiifcar la estructura global y local de batch
##################################
#Cargar paquete
devtools::install_github('theislab/kBET')

library(kBET)


# Datos sin corregir (punto de comparativa):

# Confirmar orden corecto
all(colnames(m_values_pro_684_AD) == metadata_all_filtered_357$sampleID)

# kBET - k-nearest neighbour batch effect test
batch.estimate.NO_data_corrected.batch <- kBET(t(m_values_pro_684_AD), metadata_all_filtered_357$batch)

batch.estimate.NO_data_corrected.batch$summary
#         kBET.expected kBET.observed  kBET.signif
# mean     0.03046296     0.4891667 1.165734e-16
# 2.5%     0.00000000     0.3187500 0.000000e+00
# 50%      0.02777778     0.5000000 0.000000e+00
# 97.5%    0.05555556     0.6256944 0.000000e+00

##########################################################################
# Compute a silhouette width and PCA-based measure: (Diferentes Batches)

################
## SIN CORREGIR
################
#batch = batch
batch.silhouette_no_corrected <- batch_sil(pca_no_corrected_data, metadata_all_filtered_357$batch)
# batch.silhouette_no_corrected
# [1] 0.2783481



###################
# Datos CORREGIDOS:
###################

# Confirmar orden corecto
all(metadata_all_filtered_357$sampleID == colnames(data_corrected)) # [1] TRUE

# kBET - k-nearest neighbour batch effect test
#data: a matrix (rows: cells or other observations, columns: features (genes); will be transposed if necessary)
#batch: vector or factor with batch label of each cell/observation; length has to match the size of the corresponding data dimension  
batch.estimate.data_corrected.batch <- kBET(t(data_corrected), metadata_all_filtered_357$batch)

batch.estimate.data_corrected.batch$summary
#        kBET.expected kBET.observed  kBET.signif
# mean    0.010555556    0.05222222 2.778046e-01
# 2.5%    0.000000000    0.00000000 1.227691e-10
# 50%     0.009259259    0.05555556 1.710676e-01
# 97.5%   0.027777778    0.12569444 1.000000e+00

# Compute a silhouette width and PCA-based measure:
# data: a matrix (rows: samples, columns: features (genes))
# batch: vector or factor with batch label of each cell 
batch.silhouette <- batch_sil(pca_corrected_data, metadata_all_filtered_357$batch)
# batch.silhouette
# [1] -0.01826994

# CONCLUSIONES:
# 1. No hay estructura local por batch (kBET)
# 2. No hay estructura global por batch (silhouette score)
# 3. Hubo una muy muy ligera sobre-corrección (overfitting), -0.01826994, esperable por confounding
# 4. Se removió efectivamente el efecto de lote (Confirmado con PCAs por color) :)



##############################################################
# Guardar datos preprocesados (posterior integración en SGCCA)
write.table(
  data_corrected,
  file = "Methyl_data_final_m_values.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE
)

# Confirmar lectura de datos
# data_methyl_test <- read.table(
#   "Methyl_data_final_m_values.tsv",
#   header = TRUE,
#   row.names = 1,
#   check.names = FALSE
# )

# Evaluar guardado
all.equal(data_corrected, as.matrix(data_methyl_test)) # [1] TRUE :)
 


