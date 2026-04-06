##########################################################################
#                 QA Analysis illumina 450k beadchip
##########################################################################
# 23 marzo 2026
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
  theme_bw() +
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

# Color en base a Batch
pdf("pca_m_raw_filtered_batch.pdf")
plot(
  x = pca_m_raw$x[,1],
  y = pca_m_raw$x[,2],
  col = as.factor(Metadata_filtered_723$batch)
)
dev.off()

# Color en base a sample plate
pdf("pca_m_raw_filtered_SamplePlate.pdf")
plot(
  x = pca_m_raw$x[,1],
  y = pca_m_raw$x[,2],
  col = as.factor(Metadata_filtered_723$Sample_Plate)
)
dev.off()

# Color en base a sentrixID
pdf("pca_m_raw_filtered_SentrixID.pdf")
plot(
  x = pca_m_raw$x[,1],
  y = pca_m_raw$x[,2],
  col = as.factor(Metadata_filtered_723$Sentrix_ID)
)
dev.off()

# Color en base a disease (ceradsc)
pdf("pca_m_raw_filtered_disease.pdf")
plot(
  x = pca_m_raw$x[,1],
  y = pca_m_raw$x[,2],
  col = as.factor(Metadata_filtered_723$ceradsc)
)
dev.off()

# Color en base a disease (braak)
pdf("pca_m_raw_filtered_disease_braak.pdf")
plot(
  x = pca_m_raw$x[,1],
  y = pca_m_raw$x[,2],
  col = as.factor(Metadata_filtered_723$braaksc)
)
dev.off()

# Color en base al sexo
pdf("pca_m_raw_filtered_msex.pdf")
plot(
  x = pca_m_raw$x[,1],
  y = pca_m_raw$x[,2],
  col = as.factor(Metadata_filtered_723$msex)
)
dev.off()

# Color en base a study: ROS-MAP
pdf("pca_m_raw_filtered_study.pdf")
plot(
  x = pca_m_raw$x[,1],
  y = pca_m_raw$x[,2],
  col = as.factor(Metadata_filtered_723$Study)
)
dev.off()



#
saveRDS(object = pca_m_raw, file = "pca_m_raw.rds", compress = "gzip")
saveRDS(object = idats_raw, file = "idats_raw.rds", compress = "gzip")









# Remover Batch effect (test)
library(sva)

batch <- metadata_filtered$batch


mod <- NULL

m_values_corrected <- ComBat(
  dat = m_values_filtered[, valid],   # solo las 731 muestras
  batch = batch,
  mod = mod,
  par.prior = TRUE,
  prior.plots = FALSE
)


# PCA datos corregidos
pca_noBatch <- prcomp(
                  x = t(m_values_corrected), 
             scale. = TRUE)



# Volver a graficar
pdf("pca_m_corrected.pdf")
plot(
  x = pca_noBatch$x[, 1],
  y = pca_noBatch$x[, 2]
)
dev.off()


pdf("pca_corrected_colBatch.pdf")
plot(
  x = pca_noBatch$x[, 1],
  y = pca_noBatch$x[, 2],
col = as.factor(metadata_filtered$batch)
)
dev.off()

pdf("pca_corrected_colSamplePlate.pdf")
plot(
  x = pca_noBatch$x[, 1],
  y = pca_noBatch$x[, 2],
  col = as.factor(metadata_filtered$Sample_Plate)
)
dev.off()

pdf("pca_corrected_colSentrixID.pdf")
plot(
  x = pca_noBatch$x[, 1],
  y = pca_noBatch$x[, 2],
  col = as.factor(metadata_filtered$Sentrix_ID)
)
dev.off()


#



#

































##################################
# 450k BeadChip QC & Preprocessing
##################################

#Cargar librerías
library(sesame)
library(sesameData)
library(ExperimentHub)
library(dplyr)
library(parallel)
library(BiocParallel)

# sesameDataCache("idatSignature")
# sesameDataCache()

#######################################
#           QUALITY CONTROL
#######################################

#########################
# Lectura Archivos IDAT
#########################

# # Descomentar para análisis formal !!
#
# # Extraer rutas de los archivos idat
# all_prefixes <- sesame::searchIDATprefixes(
#                     dir.name = "/STORAGE/csbig/multiomics-Arturo/methyl_data/"
#                 )
# 
# # ¿Cuántos archivos son?
# length(all_prefixes) # [1] 739 // x2 = 1478 archivos!

#Lectura de archivos IDAT


# TEST: UN SOLO CHIP
# Probar con una solo chip:
prefixes <- sesame::searchIDATprefixes(
              dir.name = "/STORAGE/csbig/multiomics-Arturo/methyl_data/5772325072/")

# Lectura de archivos IDAT correspondientes al chip 5772325072
array_5772325072 <- bplapply(
                          X = prefixes,
                        FUN = readIDATpair, 
                    BPPARAM = MulticoreParam(workers = 20)
)

# ¿Cuántos elementos tiene?
length(array_5772325072) # [1] 6 // Son 6 muestras en este beadchip


#########################
# Quality control Metrics
#########################

# 1. Bisulfite conversion
#########################

# Nota: 
# The closer the score to 1.0, the more complete the bisulfite conversion.
bis_conversion <- bplapply(
                    X = array_5772325072, 
                  FUN = sesame::bisConversionControl, 
              BPPARAM = MulticoreParam(workers = 20)
              )

# ¿La conversión de bisulfito en alguna muestra es menor que 1?
sum(!bis_conversion > 1) # [1] 0, Conversión de bisulftio bueno :)



# 2. Calcular métricas/estadísticas de control de
#     calidad para cada muestra
########################################
#
# La estadísticas incluyen:
# --Detección de sondas (probe detection)
# --Intensidades medias de señal (signal intensity)
# --Número de sondas (Number of probes)
# --Color del canal (Color Channel)
# --Sesgos de fluoróforos Cys3/Cys5 (Dye bias)
# --Métricas con valores beta (Beta value)
#
# El resultado es  una lista de métricas/estadísticas por muestra.


# 2.1 Calcula estadísticas de control de calidad para cada muestra
QC_stats_array_5772325072 <- bplapply(
                                X = array_5772325072,
                              FUN = sesame::sesameQC_calcStats,
                          BPPARAM = MulticoreParam(workers = 20)
                          )

#Observar primer elemento de la lista
QC_stats_array_5772325072[[1]] # Estadísticas de la muestra 1


# 2.2. Extraer los resultados de QC
get_stats_array_5772325072 <- lapply(
                                 X = QC_stats_array_5772325072,
                               FUN = sesameQC_getStats
                              )

get_stats_array_5772325072[1]


# 2.3 Extraer métricas de utilidad, para cada muestra

# Métricas de utilidad (selección personal)
vec1 <- c(           "num_dt",       "frac_dt", "mean_intensity",
          "mean_intensity_MU",       "mean_ii",   "mean_inb_grn", 
               "mean_inb_red",  "mean_oob_grn",   "mean_oob_red",
                 "num_probes", "num_probes_II",  "num_probes_IR", 
              "num_probes_IG",          "medR",           "medG",
                "frac_unmeth",     "frac_meth",         "num_na")

# Función para extraer métricas de utilidad
extract_metrics <- function(sample_stat){
                      sample_stat[][vec1]
                  }

beadchip_statistics <- do.call(
                          rbind,
                          bplapply(
                              X = get_stats_array_5772325072, 
                            FUN = extract_metrics, 
                        BPPARAM = MulticoreParam(workers = 20)
                          )
                       )
beadchip_statistics <- as.data.frame(beadchip_statistics)

View(beadchip_statistics)

beadchip_statistics$num_dt <- as.numeric(beadchip_statistics$num_dt)
beadchip_statistics$frac_dt <- as.numeric(beadchip_statistics$frac_dt)


# 3. Gráfico de proporción de sondas exitosas por muestra
pdf("beadchip_stats.pdf")
hist(beadchip_statistics$frac_dt)
dev.off()

# 4. Gráfico Dye bias Q-Q plot
pdf("qqplots_raw_data.pdf")
lapply(
   X = array_5772325072, 
 FUN = sesame::sesameQC_plotRedGrnQQ,
)
dev.off()

# 5.Distribución de valores beta crudos
pdf("Histogram_b_values_raw.pdf")
apply(
  X = raw_betas,
  MARGIN = 2,
  FUN = function(x){
    hist(
      x = x,
      main = "Distribución de Valores Beta",
      xlab = "Valores Beta",
      ylab = "Frecuencia")}
)
dev.off()

#AQUI TE QUEDASTEEE!!!

# 6. PCA: Identificar efectos de lote, outliers, etc

# 6.1 Extraer betas crudos
raw_betas <- do.call(
                cbind, 
                lapply(
                    X = array_5772325072, 
                  FUN = sesame::getBetas
                ) 
              )
# Extraer varianzas por sonda
var_probes <- apply(
                    X = t(raw_betas),
               MARGIN = 2,
                  FUN = var,
                na.rm = TRUE)

# Filtrar sondas con var ~ 0
probes_filtered <- t(raw_betas)[, var_probes > 0.01]

#Existen NAS, Nan, o inf ?
sum(is.na(probes_filtered)) # [1] 12
sum(is.nan(probes_filtered)) # [1] 0
sum(is.infinite(probes_filtered)) # [1] 0

# Filtrar Nas
probes_filtered <- probes_filtered[, apply(probes_filtered, 2, function(x) all(is.finite(x)))]

# PCA
pca_raw_array_5772325072 <- prcomp(
                                x = probes_filtered,
                            scale = TRUE
)

pdf("pca_raw_array_5772325072.pdf")
plot(
  x = pca_raw_array_5772325072$x[ ,1],
  y =  pca_raw_array_5772325072$x[ ,2]
)
dev.off()




# Requiere paquete "pals"
library(pals)
plot_intenseVSBetas <- lapply(
                          X = array_5772325072,
                        FUN = sesameQC_plotIntensVsBetas
                       )

sesameQC_plotIntensVsBetas(sdf = array_5772325072[[1]])



#######################################
#       Preproccesing based in QC
#######################################

# Basado en: 
# https://bioconductor.org/packages/release/bioc/vignettes/sesame/inst/doc/sesame.html#Data_Preprocessing

# 1. Quality Mask: Marcar probes problemáticos
array_5772325072_Q <- bplapply(
                            X = array_5772325072,
                          FUN = sesame::qualityMask, 
                      BPPARAM = MulticoreParam(workers = 20)
                      )

sum(array_5772325072_Q[[1]]$mask) # [1] 64144 probes masked para la primer muestra

# 2. Infer Color Channel: Inferir color del canal
array_5772325072_QC <- bplapply(
                            X = array_5772325072_Q,
                          FUN = sesame::inferInfiniumIChannel, 
                      BPPARAM = MulticoreParam(workers = 20)
                      )

# 3. Dye bias correction: Corrección de sesgo Cys3/Cys5
array_5772325072_QCD <-  bplapply(
                            X = array_5772325072_QC,
                          FUN = sesame::dyeBiasNL, 
                      BPPARAM = MulticoreParam(workers = 20)
)

# 4. pOOBAH: Detecta sondas que no supera la intensidad de background
array_5772325072_QCDP <-   bplapply(
                            X = array_5772325072_QCD,
                          FUN = sesame::pOOBAH, 
                      BPPARAM = MulticoreParam(workers = 20)
)


# 5. Noob: Background substraction (identificados en pOOBAH y QualityMask)
# Creo que tmbn normaliza (revisar)
array_5772325072_QCDPB <- bplapply(
                            X = array_5772325072_QCDP,
                          FUN = sesame::noob, 
                      BPPARAM = MulticoreParam(workers = 20)
                      )

# Extraer valores beta

processed_betas <-do.call(
                    cbind, 
                    lapply(
                        X = array_5772325072_QCDPB, 
                      FUN = sesame::getBetas
                    ) 
                  )

#######################################
#        QUALITY CONTROL after QC
#######################################

# 1. Bisulfite conversion
#########################
bis_conversion_afterQC <- bplapply(
                            X = array_5772325072_QCDPB, 
                          FUN = sesame::bisConversionControl, 
                      BPPARAM = MulticoreParam(workers = 20)
                      )

# 2.1 Calcula estadísticas de control de calidad para cada muestra
QC_stats_array_5772325072_afterQC <- bplapply(
                                        X = array_5772325072_QCDPB,
                                      FUN = sesame::sesameQC_calcStats,
                                  BPPARAM = MulticoreParam(workers = 20)
                                  )

QC_stats_array_5772325072_afterQC[1]
QC_stats_array_5772325072[1]

# 4. Gráfico Dye bias Q-Q plot
pdf("qqplots_processed_data.pdf")
lapply(
  X = array_5772325072_QCDPB, 
  FUN = sesame::sesameQC_plotRedGrnQQ,
)
dev.off()



# 5.Distriubución de valores beta
pdf("Histogram_b_values_processed.pdf")
apply(
       X = processed_betas,
  MARGIN = 2,
     FUN = function(x){
       hist(
            x = x,
         main = "Distribución de Valores Beta", 
         xlab = "Valores Beta", 
         ylab = "Frecuencia")}
)
dev.off()


# PCA data processed

# Extraer varianzas por sonda
var_probes_processed <- apply(
                          X = t(processed_betas),
                     MARGIN = 2,
                        FUN = var,
                      na.rm = TRUE
                     )

# Filtrar sondas con var ~ 0
probes_filtered <- t(processed_betas)[, var_probes_processed > 0.01]

#Existen NAS, Nan, o inf ?
sum(is.na(probes_filtered)) # [1] 12
sum(is.nan(probes_filtered)) # [1] 0
sum(is.infinite(probes_filtered)) # [1] 0

# Filtrar Nas
probes_filtered <- probes_filtered[, apply(probes_filtered, 2, function(x) all(is.finite(x)))]

# PCA
pca_processed_betas <- prcomp(
  x = probes_filtered,
  scale = TRUE
)

pdf("pca_processed_betas.pdf")
plot(
  x = pca_processed_betas$x[ ,1],
  y =  pca_processed_betas$x[ ,2]
)
dev.off()


# Generar metada_assay_bio_clinic

dim(metadata_assay) # [1] 748  12
dim(metadata_biospecimen) # [1] 748  20
dim(metadata_clinic) # [1] 3584   18











# A) mapeo metadata <--> IDAT: Se emparejan en base a la columna targetID

# 
# str(array_5772325072)
# 
# 
# prefixes[1] # "5772325072_R01C02" // Primera muestra
# 
# # Extraer el primer elemento de array // Facilitar manipulación
# sample_5772325072_R01C02 <- array_5772325072[["5772325072_R01C02"]]
# 
# # Observar data
# View(sample_5772325072_R01C02)
# 
# # Dimensiones
# dim(sample_5772325072_R01C02)
# # [1] 486427      7 
# 
# # No concuerda con los IDs de metadata
# length(metadata$TargetID) # [1] 420132
# 
# # ¿Cuáles no están en metadata y sí en idat?
# 
# metadata_ID <- metadata |>
#   select(TargetID) |>
#   unlist() |>
#   as.vector()
# 
# only_in_metadata <- sample_5772325072_R01C02 |>
#   filter(Probe_ID %in% metadata_ID)
# 
# # Estructura de la tr
# str(sample_5772325072_R01C02)
# # Classes ‘SigDF’ and 'data.frame':	486427 obs. of  7 variables:
# #   $ Probe_ID: chr  "cg00000029" "cg00000108" "cg00000109" "cg00000165" ...
# # $ MG      : int  NA NA NA NA NA NA NA NA NA 130 ...
# # $ MR      : int  NA NA NA NA NA NA NA NA NA 338 ...
# # $ UG      : int  2646 8171 3130 1200 4171 1666 6807 6369 1974 1334 ...
# # $ UR      : int  2239 504 607 3766 645 654 2757 3377 8289 13714 ...
# # $ col     : Factor w/ 3 levels "G","R","2": 3 3 3 3 3 3 3 3 3 2 ...
# # $ mask    : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
# # - attr(*, "msg")= chr "[2026-03-11 16:49:11.283059] IDAT platform: HM450"
# # - attr(*, "platform")= chr "HM450"
# # - attr(*, "controls")='data.frame':	848 obs. of  6 variables:
# #   ..$ G   : int [1:848] 100 100 12493 100 517 398 29816 30168 27986 14795 ...
# # ..$ R   : int [1:848] 19122 100 100 100 29898 31277 1818 1669 1811 606 ...
# # ..$ col : int [1:848] 64 53 31 43 41 45 35 41 45 37 ...
# # ..$ type: int [1:848] 64 53 31 43 41 45 35 41 45 37 ...
# # ..$ NA  : chr [1:848] "Red" "Purple" "Green" "Blue" ...
# # ..$ NA  : chr [1:848] "STAINING" "STAINING" "STAINING" "STAINING" ...
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
# 
# 
# sum(sample_5772325072_R01C02$mask) # [1] 0 // Se marcaron  sondas como problemáticas
# 
# unique(sample_5772325072_R01C02$col)
# 
# canal_2 <- sample_5772325072_R01C02 |>
#   filter(col == 2) |>
#   mutate(
#     M = UG, # Para este caso, probablemente así sea en ROSMAP
#     U = UR, # Para este caso, probablemente así sea en ROSMAP
#     Beta = M / (M + U)
#   )
# 
# betas_2 <- sesame::getBetas(sdf = canal_2)










