##########################################################################
#                 QA Analysis illumina 450k beadchip
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


# Flujo de trabajo (Workflow):
##############################
# 1. LECTURA DE ARCHIVOS IDAT
# 1.1 Extraer rutas de archivos (searchIDATprefixes())
# 1.2 Lectura de archivos (readIDATpair())
#
# 2. CONTROL DE CALIDAD
# 2.1 Evaluación de conversión a bisulfito (bisConversionControl())
# 2.2 Obtención de métricas de calidad (sesameQC_calcStats())
# 2.3 Análisis de datos exploratorios de las métricas de calidad (EDA)
# 2.4 Evaluación de sesgos de fluoróforos Cys3/Cys5 (sesameQC_plotRedGrnQQ())
# 2.4 Evaluación de sondas exitosas (sesameQC_plotBar())
# 2.5 Evaluación de distribución de valores beta
# 2.6 PCA de valores m
#############################################################################

#############################
# 1. LECTURA DE ARCHIVOS IDAT
#############################

##### 1.1 Extraer rutas de archivos (searchIDATprefixes())

# Extraer rutas de archivos
prefixes <- sesame::searchIDATprefixes(dir.name = "/STORAGE/csbig/multiomics-Arturo/methyl_data/")

# ¿Cuántas muestras son?
length(prefixes) # [1] 739

##### 1.2 Lectura de archivos (readIDATpair())

# Lectura de archivos idats
idats_raw <- bplapply(
                X = prefixes,
              FUN = sesame::readIDATpair,
  BPPARAM = MulticoreParam(workers = 40) # linux only
) # 34.6 GB!!

# Confirmar número de muestras
length(idats_raw) # [1] 739

# ¿Qué clase son los datos?
class(idats_raw[[1]]) # [1] "SigDF"      "data.frame"

# ¿Cómo es un SigDF?
head(idats_raw[[1]])
#     Probe_ID MG MR   UG   UR col  mask
# 1 cg00000029 NA NA 2646 2239   2 FALSE
# 2 cg00000108 NA NA 8171  504   2 FALSE
# 3 cg00000109 NA NA 3130  607   2 FALSE
# 4 cg00000165 NA NA 1200 3766   2 FALSE
# 5 cg00000236 NA NA 4171  645   2 FALSE
# 6 cg00000289 NA NA 1666  654   2 FALSE
#
# # Descripción:
# # La tabla refleja las intensidades crudas por sonda GpG. 
# # Cada fila es una sonda (probe) que mide la metilación en un sitio CpG.
# 
# # Descripción de columna:
# # Probe_ID[char]: Identificador de la sonda. Ejemplos: cg00000029, cg00000108 (cada ID corresponde a una posición genómica específica, definida en metadata)
# # MG[int]: Intensidad de fluorescencia para la sonda metilada medida en el canal verde.
# # MR[int]: Intensidad de fluorescencia para la sonda metilada medida en el canal rojo.
# # UG[int]: Intensidad de fluorescencia para la sonda no metilada medida en el canal verde.
# # UR[int]: Intensidad de fluorescencia para la sonda no metilada medida en el canal rojo.
# # col[factor]: Indica qué canal se usa para esa sonda (principal). Ejemplo: "G", "R", "2"
# # mask[logic]: Indica si la sonda fue marcada como problemática.


########################
# 2. CONTROL DE CALIDAD
########################

##### 2.1 Evaluación de conversión a bisulfito (bisConversionControl())

# lapply es necesario para aplicar a cada elemento de la lista 
bis_conversion_results <- lapply(
                            X = idats_raw,
                          FUN = bisConversionControl
)

# Graficación resultados de conversión
bis_conversion_results <- as.numeric(bis_conversion_results[1:length(bis_conversion_results)])
bis_conversion_results <- as.data.frame(bis_conversion_results)

pdf("histogram_bis_conversions.pdf")
ggplot(
  data = bis_conversion_results,
  mapping = aes(x = bis_conversion_results)
) + 
  geom_histogram(binwidth = .003) +
  labs(title = "Bisulfite Conversions") +
  xlab(label = "Bisulfite Conversions")
dev.off()

pdf("Bis_conversions_results.pdf")
ggplot(
  data = bis_conversion_results,
  mapping = aes(x = bis_conversion_results)
) + 
  geom_density() +
  scale_fill_viridis_d() +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    axis.line = element_line(color = "black")
  ) +
  labs(
    x = "Bisulfite conversion",
    y = "Density (Number of samples)",
    title = "Distribution of Bisulfite Conversion"
  ) +
  xlim(c(0.95,1.1)) +
  geom_vline(xintercept = 1,linetype = "dashed", col = "red")
dev.off()



# ¿Cuántas muestras exceden 1.1?
sum(bis_conversion_results > 1.1) # [1] 1

#¿Cuál es la muestra potencialmente a eliminar?
which(bis_conversion_results > 1.1) # [1] 453

names(idats_raw[453]) # [1] "5822038006_R03C01"



##### 2.2 Detection p-value
# Es un estimador de la calidad de la señal

pvalues <- bplapply(
              X = idats_raw, 
            FUN = function(sample){
                    pOOBAH(sample, return.pval = TRUE)
                  }, 
        BPPARAM = MulticoreParam(workers = 40)
) # 28.8 GB

# Convertir a tabla, muestras en columnas
pvalues <- do.call(
  cbind,
  pvalues)

# ¿Cuántos probes tiene al menos un pval > 0.01?
sum(rowMeans(pvalues > 0.01) != 0, na.rm = TRUE) # [1] 121303 // Se deben eliminar todos estos probes


# NAs en pvalues pueden asociarse a intensidades de señal
# no detectadas
sum(is.na(pvalues)) # [1] 1478

# Obtener pvalue promedio por muestra
pval_avgs_per_sample <- colMeans(pvalues, na.rm = TRUE)

# ¿Todas las muestras sobrepasan el umbral de selección (< 0.05)?
all(pval_avgs_per_sample < 0.05) # [1] TRUE // Todas las muestras pasan el filtro :) Indican buena calidad
# Se recomienda eliminar las muestras que no pasaron el filtro (No aplica)

# Graficar distribución de pvalues promedio por muestra, 
# Nos ayuda a visualizar las muestras que no pasaron el filtro
pdf("Histogram_pval_Avgs_per_sample.pdf")
hist(pval_avgs_per_sample, xlim = c(0.005, 0.06))
abline(v = 0.05, col = "red")
dev.off()

pval_avgs_per_sample_df <- as.data.frame(pval_avgs_per_sample)

pdf("Detection_pvalue_results.pdf")
ggplot(
  data = pval_avgs_per_sample_df,
  mapping = aes(x = pval_avgs_per_sample)
) + 
  geom_density() +
  scale_fill_viridis_d() +
  theme_classic(base_size = 14) +
  theme(
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
dev.off()

##### 2.3 Obtención de métricas de calidad (sesameQC_calcStats())

# Calcular métricas de calidad, incluyen:
# - Deteción de sondas
# - Intensidades de señal
# - Información de sondas 
# - Color del canal
# - Sesgos de fluoróforos
# - Estadísticas de valores beta

# Calcular métricas para cada muestra
qc_metrics_raw <- bplapply(
                      X = idats_raw,
                    FUN = sesame::sesameQC_calcStats,
                BPPARAM = MulticoreParam(workers = 40)
) # 16.7 MB

# Exraer métricas/ estadísticas
stats <- lapply(
            X = qc_metrics_raw,
          FUN = sesameQC_getStats
) # 6.7 MB

# Extraer nombres de las métricas
metric_names <- names(stats[[1]])

# Extraer métricas para cada muestra y convertirlo a una tabla
qc_statistics <- do.call(
                    rbind,
                    lapply(
                        X = stats, 
                      FUN = function(sample_stat){
                              sample_stat[][metric_names]
                            }
                    )
                  )

# Convertir tabla a data frame
qc_statistics <- as.data.frame(qc_statistics)

View(qc_statistics)



##### 2.4 Evaluación de sesgos de fluoróforos Cys3/Cys5 (sesameQC_plotRedGrnQQ())

# Generar gráficas qqplot de las intensidades de señal roja y verde
pdf("raw_qqplots_intensityRedGrn.pdf")
lapply(
  X = idats_raw, 
  FUN = sesameQC_plotRedGrnQQ,
)
dev.off()

# Resumen qqplots
pdf("RGratio_vs_RGdistortion_raw.pdf")
plot(qc_statistics$RGratio, qc_statistics$RGdistort,
     xlab="RG ratio (shift)",
     ylab="RG distortion (curvature)",
     pch=16)
abline(h=0, lty=2)
abline(v=1, lty=2)
dev.off()


##### 2.5 Evaluación de sondas exitosas (sesameQC_plotBar())
pdf("raw_plotBar.pdf")
sesameQC_plotBar(
  lapply(
    X = idats_raw,
    FUN = sesameQC_calcStats,
    "detection"
  )
)
dev.off()


qc_statistics$frac_dt <- as.numeric(qc_statistics$frac_dt)
succes_probes <- qc_statistics %>%
                  select(frac_dt) %>%
                  mutate(frac_dt_perc = frac_dt * 100)



pdf("Probe_Detection_Succes_results.pdf")
ggplot(
  data = succes_probes,
  mapping = aes(x = frac_dt_perc)
) + 
  geom_density(aes(y = after_stat(count))) +
  scale_fill_viridis_d() +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    axis.line = element_line(color = "black")
  ) +
  labs(
    x = "Detection Succes (%)",
    y = "Density (Number of samples)",
    title = "Probe Detection Succes"
  ) +
  xlim(c(0,110)) +
  geom_vline(xintercept = 95, linetype = "dashed", col = "red")
dev.off()


##### 2.6 Evaluación de distribución de valores beta

# Extraer betas crudos y transformarlo a una tabla
betas_raw <- do.call(
  cbind, 
  lapply(
    X = idats_raw, 
    FUN = function(prefix){
      sesame::getBetas(prefix)
    }
  ) 
) # 2.9 GB

# Gráfico de distribución de valores beta promedio (todas las muestras)
pdf("Dist_b_values_raw.pdf")
hist(betas_raw)
dev.off()



betas_long <- betas_raw %>%
  as.data.frame() %>%
  rownames_to_column(var = "CpG") %>%
  pivot_longer(
    cols = -CpG,
    names_to = "sample",
    values_to = "beta"
  )

pdf("b_values_raw_results.pdf")
ggplot(betas_long, aes(x = beta, group = sample)) +
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


betas_raw_summary <- betas_long %>%
  group_by(sample) %>%
  summarise(
    mean_beta = mean(beta, na.rm = TRUE),
    median_beta = median(beta, na.rm = TRUE),
    sd_beta = sd(beta, na.rm = TRUE)
  )

pdf("scatterplot_betas_raw.pdf")
ggplot(betas_raw_summary, aes(x = mean_beta, y = sd_beta, label = sample)) +
  geom_point() +
  geom_text(size = 2, vjust = -0.5) +
  theme_minimal()
dev.off()


# Gráficos de distribución de valores beta, para cada muestra
pdf("Dist_b_values_raw_persample.pdf")
apply(
  X = betas_raw,
  MARGIN = 2,
  FUN = function(x){
    hist(
      x = x,
      main = "Distribución de Valores Beta",
      xlab = "Valores Beta",
      ylab = "Frecuencia")}
)
dev.off()


##### 2.8 Identificación de sondas a filtrar
# Sondas con NA,
# Sondas con mapeo a cromosomas sexuales, 
# Sondas con SNPs en el sitio CpG
# Sondas con hibridación no específica

## NAs
# ¿Existen sondas con NA?
sum(is.na(betas_raw)) # [1] 1478 // Provienen de la fórmula del cálculo de valores beta. fallos de detección

#¿Qué proporción de NAs hay?
mean(is.na(betas_raw)) # [1] 4.111614e-06
# Esta proporción es mínima, practicamente 0. En revisiones, se recomienda eliminarlas.
# Por tanto, con este filtrado no se necesitará imputación de datos (no es necesario, para este caso)




#### Dendograma (clustering)

euc_distbetas_raw <- hclust(dist(t(betas_raw)))

pdf("dendogram_betas_raw.pdf")
plot(euc_distbetas_raw)
dev.off()


##### 2.7 PCA de valores m

## Filtrar muestras con mapeo a metadata (necesario para agg color)

# SampleID de idats_raw
sampleID_idat_data <- names(idats_raw)

length(sampleID_idat_data) # [1] 739

# SampleID de metadata_all
sampleID_metadata_all <- metadata_all$sampleID

length(sampleID_metadata_all) # [1] 740

samples_keep <- intersect(x = sampleID_idat_data, y = sampleID_metadata_all) 

# ¿Cuántas muestras quedan?
length(samples_keep) # [1] 723 // Algunas muestras no mapean a metadata :(

# Filtrar muestras en betas_raw
betas_raw_723 <- betas_raw[,samples_keep]

dim(betas_raw_723) # [1] 486427    723



## Filtrado necesario en PCA (Eliminar sondas con NA)

# Calcular proporción no NAs por sonda
probes_keep <- rowMeans(!is.na(betas_raw_723))

# remove any probes that have failed in one or more samples
betas_raw_filtered <- betas_raw_723[probes_keep == 1, ]

dim(betas_raw_filtered) # [1] 486425    723// Se eliminaron solo 2 probes

# Confirmar que no ocurra error en PCA
sum(is.na(betas_raw_filtered)) # [1] 0
sum(is.nan(betas_raw_filtered)) # [1] 0
sum(is.infinite(betas_raw_filtered)) # [1] 0


# Convertir valores beta a valores M (mayor homocedasticidad, ~ normal)
m_values_raw <- BetaValueToMValue(
  b = betas_raw_filtered
) # 2.8 GB

# Gráficos de distribución de valores m promedio (todas las muestras)
pdf("Hist_m_values_raw_all.pdf")
hist(m_values_raw)
dev.off()


# ¿Existen NA, Nan o inf? (m-values)
sum(is.na(m_values_raw)) # [1] 0
sum(is.nan(m_values_raw)) # [1] 0
sum(is.infinite(m_values_raw)) # [1] 592 // Provienen de la fórmula del cálculo de valores m

# Extraer sondas de interés (no NA, no inf)
probes_filter <- apply(
  X = m_values_raw,
  MARGIN =  1,
  FUN =  function(x){
    all(is.finite(x))
  } 
)

length(probes_filter) # [1] 486425


# Filtrar sondas de interés
m_values_filtered <- m_values_raw[probes_filter, ]

dim(m_values_filtered) # [1] 485979    723 // En total se eliminaron 448 probes

# Realizar PCA con valores m filtrados
pca_m_raw <- prcomp(
  x = t(m_values_filtered), 
  scale. = TRUE) # 2.9 GB


##### 2.8 Graficar PCA valores m crudos
pdf("pca_m_raw_723_sampleID.pdf")
plot(
  x = pca_m_raw$x[ ,1], 
  y = pca_m_raw$x[ ,2], 
  col = as.factor(Metadata_filtered_723$sampleID)
)
dev.off()

pca_var <- pca_m_raw$sdev^2

pca_var_per <- round(pca_var / sum(pca_var) * 100, 1)

pdf("Scree_plot_pca_m_raw.pdf")
barplot(height = pca_var_per, main = "Scree plot", xlab = "Principal Component", ylab = "Percent Variation")
dev.off()

#Formatear PCA  a dataframe para ggplot
pca_m_raw_df <- data.frame(
  sample = rownames(pca_m_raw$x),
  X = pca_m_raw$x[,1],
  Y = pca_m_raw$x[,2]
)

# Graficar
pdf("PCA_raw.pdf")
ggplot(
  data = pca_m_raw_df,
  mapping = aes(x = X, y = Y, label = sample)
) +
  geom_text() +
  xlab(paste("PC1 - ", pca_var_per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca_var_per[2], "%", sep = "")) +
  theme_classic() +
  ggtitle("PCA")
dev.off()
  

# Colorear en base a distintas variables

#Filtrar metadata

sampleID_pca_m_raw <- rownames(pca_m_raw$x)

Metadata_filtered_723 <- metadata_all %>%
                          filter(sampleID %in% sampleID_pca_m_raw)

# ¿El orden de las muestras es el mismo que en metadata?
all(sampleID_pca_m_raw == Metadata_filtered_723$sampleID) # [1] FALSE // Se requiere ordenar metadata

# Ordenar metadata en base al orden de las muestras en el PCA
Metadata_filtered_723 <- Metadata_filtered_723[match(x = sampleID_pca_m_raw, table = Metadata_filtered_723$sampleID), ]

# Confirmar orden
all(sampleID_pca_m_raw == Metadata_filtered_723$sampleID) # [1] TRUE :)

all(pca_m_raw_df$sample == Metadata_filtered_723$sampleID) # [1] TRUE :)

# Color en base a Batch
pdf("pca_m_raw_filtered_batch.pdf")
pca_m_raw_df %>% 
      ggplot(mapping = aes(x = X, y = Y)
      ) +
      geom_point() +
      aes(colour = as.factor(Metadata_filtered_723$batch)) +
      scale_color_discrete(name = "Batch") +
      xlab(paste("PC1 - ", pca_var_per[1], "%", sep = "")) +
      ylab(paste("PC2 - ", pca_var_per[2], "%", sep = "")) +
      theme_classic() +
      ggtitle("PCA") +
      stat_ellipse(geom = "polygon", aes(fill = as.factor(Metadata_filtered_723$batch)), alpha = 0.2, show.legend = FALSE)
dev.off()

# Color en base a sample plate
pdf("pca_m_raw_filtered_SamplePlate.pdf")
pca_m_raw_df %>% 
  ggplot(mapping = aes(x = X, y = Y)
  ) +
  geom_point() +
  aes(colour = as.factor(Metadata_filtered_723$Sample_Plate)) +
  scale_color_discrete(name = "Sample Plate") +
  xlab(paste("PC1 - ", pca_var_per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca_var_per[2], "%", sep = "")) +
  theme_classic() +
  ggtitle("PCA") +
  stat_ellipse(geom = "polygon", aes(fill = as.factor(Metadata_filtered_723$Sample_Plate)), alpha = 0.2, show.legend = FALSE)
dev.off()

# Color en base a sentrixID
pdf("pca_m_raw_filtered_SentrixID.pdf")
pca_m_raw_df %>% 
  ggplot(mapping = aes(x = X, y = Y)
  ) +
  geom_point() +
  aes(colour = as.factor(Metadata_filtered_723$Sentrix_ID)) +
  scale_color_discrete(guide = "none") +
  xlab(paste("PC1 - ", pca_var_per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca_var_per[2], "%", sep = "")) +
  theme_classic() +
  ggtitle("PCA") +
  stat_ellipse(geom = "polygon", aes(fill = as.factor(Metadata_filtered_723$Sentrix_ID)), alpha = 0.2, show.legend = FALSE)
dev.off()

# Color en base a disease (ceradsc)
pdf("pca_m_raw_filtered_disease.pdf")
pca_m_raw_df %>% 
  ggplot(mapping = aes(x = X, y = Y)
  ) +
  geom_point() +
  aes(colour = as.factor(Metadata_filtered_723$ceradsc)) +
  scale_color_discrete(name = "Cerad") +
  xlab(paste("PC1 - ", pca_var_per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca_var_per[2], "%", sep = "")) +
  theme_classic() +
  ggtitle("PCA") +
  stat_ellipse(geom = "polygon", aes(fill = as.factor(Metadata_filtered_723$ceradsc)), alpha = 0.2, show.legend = FALSE)
dev.off()

# Color en base a disease (braak)
pdf("pca_m_raw_filtered_disease_braak.pdf")
pca_m_raw_df %>% 
  ggplot(mapping = aes(x = X, y = Y)
  ) +
  geom_point() +
  aes(colour = as.factor(Metadata_filtered_723$braaksc)) +
  scale_color_discrete(name = "Braak") +
  xlab(paste("PC1 - ", pca_var_per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca_var_per[2], "%", sep = "")) +
  theme_classic() +
  ggtitle("PCA") +
  stat_ellipse(geom = "polygon", aes(fill = as.factor(Metadata_filtered_723$braaksc)), alpha = 0.2, show.legend = FALSE)
dev.off()

# Color en base al sexo
pdf("pca_m_raw_filtered_msex.pdf")
pca_m_raw_df %>% 
  ggplot(mapping = aes(x = X, y = Y)
  ) +
  geom_point() +
  aes(colour = as.factor(Metadata_filtered_723$msex)) +
  scale_color_discrete(name = "Sex") +
  xlab(paste("PC1 - ", pca_var_per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca_var_per[2], "%", sep = "")) +
  theme_classic() +
  ggtitle("PCA") +
  stat_ellipse(geom = "polygon", aes(fill = as.factor(Metadata_filtered_723$msex)), alpha = 0.2, show.legend = FALSE)
dev.off()

# Color en base a study: ROS-MAP
pdf("pca_m_raw_filtered_study.pdf")
pca_m_raw_df %>% 
  ggplot(mapping = aes(x = X, y = Y)
  ) +
  geom_point() +
  aes(colour = as.factor(Metadata_filtered_723$Study)) +
  scale_color_discrete(name = "Study") +
  xlab(paste("PC1 - ", pca_var_per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca_var_per[2], "%", sep = "")) +
  theme_classic() +
  ggtitle("PCA") +
  stat_ellipse(geom = "polygon", aes(fill = as.factor(Metadata_filtered_723$Study)), alpha = 0.2, show.legend = FALSE)
dev.off()


# Fin QA
