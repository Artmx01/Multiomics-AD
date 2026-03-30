##########################################################################
#                 Quality Control illumina 450k beadchip
##########################################################################
# 23 marzo 2026
# Arturo, BM

# Cargar librerias
library(sesame)
library(sesameData)
library(BiocParallel)
library(dplyr)
library(ggplot2)
library(illuminaio)


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

# En la gráfica se observa una muestra mayor que 1.1. Esto es indicativo de una posible
# mala calidad de muestra. Es necesario confirmar con detection p-value 

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
sum(is.na(betas_raw)) # [1] 1478 // Provienen de la fórmula del cálculo de valores beta

#¿Qué proporción de NAs hay?
mean(is.na(betas_raw)) # [1] 4.111614e-06
# Esta proporción es mínima, practicamente 0. En revisiones, se recomienda eliminarlas.
# Por tanto, con este filtrado no se necesitará imputación de datos (no es necesario, para este caso)


## Sex probes
# En el metadata no hay sondas con mapeo a cromosomas sexuales
metadata <- vroom::vroom("/STORAGE/csbig/multiomics-Arturo/methyl_data/metadata/ROSMAP_arrayMethylation_metaData.tsv")

# Sondas presentes en metadata
unique(metadata$CHR) # [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22

# Extraer ID probes de metadata
nosex_probes <- metadata %>%
  select(TargetID) %>%
  unlist() %>% 
  as.vector()

## SNPs probes
manifest <- illuminaio::readBPM(file = "/STORAGE/csbig/multiomics-Arturo/methyl_data/HumanMethylation450K/HumanMethylation450_15017482_v.1.1.bpm")


## Cross reactive probes



##### 2.7 PCA de valores m

## Filtrado necesario en PCA (Eliminar sondas con NA)

# Calcular proporción no NAs por sonda
keep <- rowMeans(!is.na(betas_raw))

# remove any probes that have failed in one or more samples
betas_raw_filtered <- betas_raw[keep == 1, ] 

dim(betas_raw_filtered) # [1] 486425    739 // Se eliminaron solo 2 probes

# Confirmar que no ocurra error en PCA
sum(is.na(betas_raw_filtered)) # [1] 0
sum(is.nan(betas_raw_filtered)) # [1] 0
sum(is.infinite(betas_raw_filtered)) # [1] 0


# Convertir valores beta a valores M (mayor homocedasticidad, ~ normal)
m_values_raw <- BetaValueToMValue(
  b = betas_raw_filtered
) # 2.9 GB

# Gráficos de distribución de valores m promedio (todas las muestras)
pdf("Hist_m_values_raw_all.pdf")
hist(m_values_raw)
dev.off()


# ¿Existen NA, Nan o inf? (m-values)
sum(is.na(m_values_raw)) # [1] 0
sum(is.nan(m_values_raw)) # [1] 0
sum(is.infinite(m_values_raw)) # [1] 603 // Provienen de la fórmula del cálculo de valores m

# Extraer sondas de interés (no NA, no inf)
probes_filter <- apply(
  X = m_values_raw,
  MARGIN =  1,
  FUN =  function(x){
    all(is.finite(x))
  } 
)

length(probes_filter)


# Filtrar sondas de interés
m_values_filtered <- m_values_raw[probes_filter, ]

dim(m_values_filtered) # [1] 485971    739 // En total se eliminaron 454 probes

# Realizar PCA con valores m filtrados
pca_m_raw <- prcomp(
  x = t(m_values_filtered), 
  scale. = TRUE)


##### 2.8 Graficar PCA valores m crudos
pdf("pca_m_raw.pdf")
plot(
  x = pca_m_raw$x[ ,1], 
  y = pca_m_raw$x[ ,2]
)
dev.off()



metadata_assay <- vroom::vroom(file = "/STORAGE/csbig/multiomics-Arturo/methyl_data/metadata/ROSMAP_assay_methylationArray_metadata.csv")
View(metadata_assay)


samples_order <- rownames(pca_m_raw$x)

metadata_assay <- metadata_assay %>%
  mutate(
    sampleID = paste(Sentrix_ID, Sentrix_Row_Column, sep = "_")
  )

idx <- match(samples_order, metadata_assay$sampleID)
valid <- !is.na(idx)

pca_filtered <- pca_m_raw$x[valid, ]
metadata_filtered <- metadata_assay[idx[valid], ]

all(rownames(pca_filtered) == metadata_filtered$sampleID)

pdf("pca_m_raw_filtered.pdf")
plot(
  x = pca_filtered[,1],
  y = pca_filtered[,2],
  col = as.factor(metadata_filtered$batch)
)
dev.off()

# Color en base a Batch
pdf("pca_m_raw_colBatch.pdf")
plot(
  x = pca_filtered[,1],
  y = pca_filtered[,2],
  col = as.factor(metadata_filtered$batch)
)
dev.off()

# Color en base a sample plate
pdf("pca_m_raw_colSamplePlate.pdf")
plot(
  x = pca_filtered[,1],
  y = pca_filtered[,2],
  col = as.factor(metadata_filtered$Sample_Plate)
)
dev.off()

# Color en base a sentrixID
pdf("pca_m_raw_colSentrixID.pdf")
plot(
  x = pca_filtered[,1],
  y = pca_filtered[,2],
  col = as.factor(metadata_filtered$Sentrix_ID)
)
dev.off()