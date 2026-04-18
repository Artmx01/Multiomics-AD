##########################################################################
#                 QC Analysis illumina 450k beadchip
##########################################################################
# 17 abril 2026
# Arturo, BM


# Cargar librerias
library(sesame)
library(sesameData)
library(BiocParallel)
library(ggplot2)
library(illuminaio)
library(tidyverse)
library(sva)

# Previamente, revisar methyl450k_QA.R


##########################################
# FILTRADO DE MUESTRAS
##########################################

## 1. Filtrar muestras solo con mapeo a metadata

# SampleID de idats_raw
sampleID_idat_data <- names(idats_raw)

length(sampleID_idat_data) # [1] 739

# SampleID de metadata_all
sampleID_metadata_all <- metadata_all$sampleID

length(sampleID_metadata_all) # [1] 740

samples_keep <- intersect(x = sampleID_idat_data, y = sampleID_metadata_all) 

# ¿Cuántas muestras quedan?
length(samples_keep) # [1] 723 // Algunas muestras no mapean a metadata :(

# Filtrar muestras en idats_raw
idats_raw_723 <- idats_raw[samples_keep]

# Confirmar filtrado
length(idats_raw_723) # [1] 723


## 2. Filtrar muestras con conversión a bisulfito < 0.7
# Todas las muestras pasaron el filtro :)

## 3. Filtrar muestras con Detection p-value > 0.05
# Todas las muestras pasaron el filtro :)

## 4. Filtrar muestras con Succes probes Detection < 95%
qc_statistics_filtered <- qc_statistics[samples_keep, ]

# Columna frac_dt contiene la proporción de sondas exitosas
samples_keep <- qc_statistics_filtered %>%
                  filter(frac_dt >= 0.95) %>%
                  rownames()

length(samples_keep) # [1] 688

# Filtrar muestras
idats_raw_688 <- idats_raw_723[samples_keep]

# Confirmar filtrado
length(idats_raw_688) # [1] 688 // Se eliminan 35 muestras

## 5. Filtrado de muestras atípicas (outliers)

# Muestras identificadas con plot promedio_de_betas vs desvest o PCA
any(names(idats_raw_688) == "5815381015_R06C02") # [1] TRUE
any(names(idats_raw_688) == "5822038005_R03C02") # [1] FALSE
any(names(idats_raw_688) == "5815381027_R05C01") # [1] FALSE
any(names(idats_raw_688) == "5822038006_R05C01") # [1] FALSE
any(names(idats_raw_688) == "5822071004_R01C02") # [1] TRUE
any(names(idats_raw_688) == "5822020001_R05C01") # [1] TRUE

no_keep <- c("5815381015_R06C02", "5822071004_R01C02", "5822020001_R05C01") # Solo TRUE

samples_688 <- names(idats_raw_688)

samples_keep <- samples_688[!samples_688 %in% no_keep]

# Filtrar muestras
idats_raw_685 <- idats_raw_688[samples_keep]

# Confirmar filtrado
length(idats_raw_685) # [1] 685 // Se eliminan 3 muestras


##########################################
# SESAME PREPROCESSING
##########################################

## Detectar/Enmascarar probes "ruidosos" (QualityMask())
## Inferir color del canal (inferInfiniumIChannel())
## Corregir sesgos de fluoróforos (Cys3/Cys5) (DybiasNL())
## Identificar probes fallidos pOOBAH (pOOBAH())
## Reducción de ruido de fondo/ Background substraction (noob())

idats_pro_685 <- bplapply(
                    X = idats_raw_685, 
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

# Extraer betas y transformarlo a matrix
betas_pro_685 <- do.call(
                  cbind, 
                  lapply(
                      X = idats_pro_685, 
                    FUN = function(prefix){
                            sesame::getBetas(prefix)
                          }
                   ) 
                 )


##########################################
# FILTRADO DE SONDAS
##########################################

## 1. Sondas con valores faltantes (NA, NaN, infinite)

# ¿Existen valores faltantes?
sum(is.na(betas_pro_685)) # [1] 50619771
sum(is.nan(betas_pro_685)) # [1] 0
sum(is.infinite(betas_pro_685)) # [1] 0

# ¿Cual es la proporción de NAs?
mean(is.na(betas_pro_685)) # [1] 0.1519189 // 15% de los datos

# ¿Cuántas sondas tienen datos completos (no NAs)?
sum(rowMeans(!is.na(betas_pro_685)) == 1) # [1] 375607 // de 486427

# remove any probes that have failed in one or more samples
keep <- rowMeans(!is.na(betas_pro_685))

betas_pro_685_filtered <- betas_pro_685[keep == 1, ]

# Confirmar filtrado
dim(betas_pro_685_filtered) # [1] 375607    685 // se eliminaron 110820 sondas
sum(is.na(betas_pro_685_filtered)) # [1] 0
sum(is.nan(betas_pro_685_filtered)) # [1] 0
sum(is.infinite(betas_pro_685_filtered)) # [1] 0

# Identificar sondas que fallan en al menos 1% de las muestras (consideranddo pval > 0.01)
bad_probes <- rownames(pvalues)[rowMeans(pvalues > 0.01) > 0.01]

# ¿Cuántos bad probes son?
length(bad_probes) # [1] 86773 // de 375607

probes_all <- rownames(betas_pro_685_filtered)
no_keep <- intersect(x = probes_all, y = bad_probes) # [1] 27170
keep_probes <- !(probes_all %in% bad_probes)

betas_pro_685_filtered_pval <- betas_pro_685_filtered[keep_probes, ]

# Confirmar filtrado
dim(betas_pro_685_filtered_pval) # [1] 348437    685 // Se eliminan 27170 probes


## 2. Sondas con mapeo a cromosomas sexuales

# Utilizaremos el archivo metadata que contiene solo sondas con mapeo a cromosomas
# somáticos

somatic_probes <- metadata %>%
  select(TargetID) %>% 
  unlist() %>%
  as.vector()

probes_betas_pro <- rownames(betas_pro_685_filtered_pval)

keep_probes <- intersect(x = somatic_probes, y = probes_betas_pro)

# Filtrar sondas con mapeo a cromosomas sexuales
betas_pro_685_filtered_pval_nosex <- betas_pro_685_filtered_pval[keep_probes, ]

# Confirmar filtrado
dim(betas_pro_685_filtered_pval_nosex) # [1] 335819    685 // Se eliminan 12618


## 3. Sondas con SNPs en el sitio CpG, cross-reactive:

# Zhou et al, 2016
# https://pubmed.ncbi.nlm.nih.gov/27924034/
# DOI: 10.1093/nar/gkw967

url <- "https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno/HM450/archive/202209/HM450.hg19.manifest.tsv.gz"
file <- "HM450.hg19.manifest.tsv"

download.file(url = url, destfile = file)

HM450.hg19.manifest.tsv <- vroom(file = "HM450.hg19.manifest.tsv")

View(HM450.hg19.manifest.tsv)

bad_probes <- HM450.hg19.manifest.tsv %>%
  filter(MASK_general == TRUE) %>%
  select(probeID) %>%
  unlist() %>%
  as.vector()

# ¿Cuántas bad_probes son?
length(bad_probes) # [1] 60466

probes_betas_pro <- rownames(betas_pro_685_filtered_pval_nosex)

keep_probes <- !(probes_betas_pro %in% bad_probes)

betas_pro_685_filtered_pval_nosex_noSNP_noCrossreactive <- betas_pro_685_filtered_pval_nosex[keep_probes, ]

# Confirmar filtrado
dim(betas_pro_685_filtered_pval_nosex_noSNP_noCrossreactive) # [1] 335299    685 // se eliminaron 520 sondas



##########################################
# REMOVER EFECTO DE LOTE
##########################################

# Convertir beta values a m values
m_values_pro<- BetaValueToMValue(
  b = betas_pro_685_filtered_pval_nosex_noSNP_noCrossreactive
)

# Eliminar muestra 5822038012_R02C01, contiene NA en metadata
no_keep <- "5822038012_R02C01"

m_values_pro_684 <- m_values_pro[ ,!(colnames(m_values_pro) %in% no_keep)]


# También es necesario filtrar metadata
samples_684 <- colnames(m_values_pro_684)

metadata_filtered_684 <- Metadata_filtered_723 %>%
                            filter(sampleID %in% samples_684)

# Confirmar orden entre muestras y metadata
all(samples_684 == metadata_filtered_684$sampleID) # [1] TRUE :)

# Describir variables al modelo
batch <- metadata_filtered_684$batch

mod <- model.matrix(~ msex + educ + race + apoe_genotype + ceradsc, data = metadata_filtered_684)

# Remover efecto de lote
m_values_pro_684_noBatch <- ComBat(
                          dat = m_values_pro_684, 
                        batch = batch, 
                          mod = mod)




###############
# FIN :)
###############

# BALANCE DE VARIABLES

# Evaluar si cerad y batch son independientes
chisq.test(table(metadata_filtered_684$ceradsc,metadata_filtered_684$batch))
# Pearson's Chi-squared test
# 
# data:  table(metadata_filtered_684$ceradsc, metadata_filtered_684$batch)
# X-squared = 13.019, df = 3, p-value =
# 0.004596

# Conclusiones:
# No hay independencia entre cerad y batch = Existe asociación entre las variables
# Esto podría afectar a combat de la siguiente manera:
# 1. Eliminar parte de la señal biológica real
# 2. Sobre-ajustar

# Ver proporciones de grupo en cada batch
table(metadata_filtered_684$ceradsc, metadata_filtered_684$batch)
#     0   1
# 1  75 139
# 2  81 145
# 3  24  48
# 4  87  85

# Visualizar proporciones
pdf("prop_cerad_batch.pdf")
ggplot(metadata_filtered_684, aes(x = batch, fill = as.factor(ceradsc))) +
  geom_bar() +
  theme_classic()
dev.off()


table(metadata_filtered_684$batch)
# 0   1 
# 267 417 // Diferencia de 150 individuos

# Distribución del fenotipo cerad en cada lote (0 y 1):

# Batch 0
metadata_filtered_684 %>% 
  filter(batch == 0) %>% 
  select(ceradsc) %>% 
  table()
# ceradsc
# 1  2  3  4 
# 75 81 24 87 

# Batch 1
metadata_filtered_684 %>% 
  filter(batch == 1) %>% 
  select(ceradsc) %>% 
  table()
# ceradsc
# 1   2   3   4 
# 139 145  48  85 

# Distribución del sexo en cada lote (0 y 1):

# Batch 0
metadata_filtered_684 %>% 
  filter(batch == 0) %>% 
  select(msex) %>% 
  table()
# msex
# 0   1 
# 168  99 

# Batch 1
metadata_filtered_684 %>% 
  filter(batch == 1) %>% 
  select(msex) %>% 
  table()
# msex
# 0   1 
# 264 153 

# Hay más individuos en el lote 1 que en el lote 2 (150+)

# Confirmar orden entre muestras y metadata
all(metadata_filtered_684$sampleID == colnames(m_values_pro_684)) # [1] TRUE



# MODELOS DE COMBAT
#######################
# 1. Modelo nulo

# Describir variables al modelo
batch <- metadata_filtered_684$batch

mod <- NULL

# Remover efecto de lote, modelo nulo
test_combat_null <- ComBat(
  dat = m_values_pro_684, 
  batch = batch, 
  mod = mod)

# PCA
pca_null <- prcomp(
  x = t(test_combat_null), 
  scale. = TRUE)

# Graficar
all(rownames(pca_null$x) == metadata_filtered_684$sampleID)

pdf("pca_null.pdf")
plot(
  x = pca_null$x[, 1],
  y = pca_null$x[, 2],
  col = as.factor(metadata_filtered_684$batch)
)
dev.off()

#######################
# 2. Modelo: ~ msex

# Describir variables al modelo
batch <- metadata_filtered_684$batch

mod <- model.matrix(~ msex, metadata_filtered_684)

# Remover efecto de lote, modelo ~ msex
test_combat_msex <- ComBat(
  dat = m_values_pro_684, 
  batch = batch, 
  mod = mod)

# PCA
pca_msex <- prcomp(
  x = t(test_combat_msex), 
  scale. = TRUE)

# Graficar
all(rownames(pca_msex$x) == metadata_filtered_684$sampleID)

pdf("pca_msex.pdf")
plot(
  x = pca_msex$x[, 1],
  y = pca_msex$x[, 2],
  col = as.factor(metadata_filtered_684$batch)
)
dev.off()

############################
# 3. Modelo: ~ educ

# Describir variables al modelo
batch <- metadata_filtered_684$batch

mod <- model.matrix(~ educ, metadata_filtered_684)

# Remover efecto de lote, modelo ~ educ
test_combat_educ <- ComBat(
  dat = m_values_pro_684, 
  batch = batch, 
  mod = mod)

# PCA
pca_educ <- prcomp(
  x = t(test_combat_educ), 
  scale. = TRUE)

# Graficar
all(rownames(pca_educ$x) == metadata_filtered_684$sampleID)

pdf("pca_educ.pdf")
plot(
  x = pca_educ$x[, 1],
  y = pca_educ$x[, 2],
  col = as.factor(metadata_filtered_684$batch)
)
dev.off()


############################
# 3. Modelo: ~ ceradsc

# Describir variables al modelo
batch <- metadata_filtered_684$batch

mod <- model.matrix(~ as.factor(ceradsc), metadata_filtered_684)

# Remover efecto de lote, modelo ~ educ
test_combat_cerad <- ComBat(
  dat = m_values_pro_684, 
  batch = batch, 
  mod = mod)

# PCA
pca_cerad <- prcomp(
  x = t(test_combat_cerad), 
  scale. = TRUE)

# Graficar
all(rownames(pca_cerad$x) == metadata_filtered_684$sampleID)

pdf("pca_cerad.pdf")
plot(
  x = pca_cerad$x[, 1],
  y = pca_cerad$x[, 2],
  col = as.factor(metadata_filtered_684$batch)
)
dev.off()


##########################
# 4. Modelo: ~ msex + educ

# Describir variables al modelo
batch <- metadata_filtered_684$batch

mod <- model.matrix(~ msex + educ, metadata_filtered_684)

# Remover efecto de lote, modelo ~ msex + educ
test_combat_msex_educ <- ComBat(
  dat = m_values_pro_684, 
  batch = batch, 
  mod = mod)

# PCA
pca_msex_educ <- prcomp(
  x = t(test_combat_msex_educ), 
  scale. = TRUE)

# Graficar
all(rownames(pca_msex_educ$x) == metadata_filtered_684$sampleID)

pdf("pca_msex_educ.pdf")
plot(
  x = pca_msex_educ$x[, 1],
  y = pca_msex_educ$x[, 2],
  col = as.factor(metadata_filtered_684$batch)
)
dev.off()


###################
###################

batch <- interaction(metadata_filtered_684$batch, metadata_filtered_684$Sample_Plate, metadata_filtered_684$Sentrix_ID)

mod <- model.matrix(~ ceradsc, metadata_filtered_684)

# Remover efecto de lote, modelo ~ msex + educ
test_combat_no_batch_samplePlate_SentrixID <- ComBat(
  dat = m_values_pro_684, 
  batch = batch, 
  mod = mod)
# Using the 'mean only' version of ComBat
# Found1280batches
# Note: one batch has only one sample, setting mean.only=TRUE
# Adjusting for1covariate(s) or covariate level(s)
# Error en ComBat(dat = m_values_pro_684, batch = batch, mod = mod): 
#   The covariate is confounded with batch! Remove the covariate and rerun ComBat

# MODELADO DE BATCH COMO SENTRIX_ID
batch <- metadata_filtered_684$Sentrix_ID

mod <- model.matrix(~ ceradsc, metadata_filtered_684)

# Remover efecto de lote (batch = SentrixID)
test_combat_no_SentrixID <- ComBat(
  dat = m_values_pro_684, 
  batch = batch, 
  mod = mod)
# Using the 'mean only' version of ComBat
# Found64batches
# Note: one batch has only one sample, setting mean.only=TRUE
# Adjusting for1covariate(s) or covariate level(s)
# Standardizing Data across genes
# Fitting L/S model and finding priors
# Finding parametric adjustments
# Adjusting the Data


# PCA
pca_noSentrixID <- prcomp(
  x = t(test_combat_no_SentrixID), 
  scale. = TRUE)

# Graficar
all(rownames(pca_noSentrixID$x) == metadata_filtered_684$sampleID)

pdf("pca_no_SentrixID.pdf")
plot(
  x = pca_noSentrixID$x[, 1],
  y = pca_noSentrixID$x[, 2],
  col = as.factor(metadata_filtered_684$braaksc)
)
dev.off()

################################################
# MODELADO DE BATCH COMO SAMPLE_PLATE
batch <- metadata_filtered_684$Sample_Plate

mod <- model.matrix(~ ceradsc, metadata_filtered_684)

# Remover efecto de lote (batch = SamplePlate)
test_combat_no_SamplePlate <- ComBat(
  dat = m_values_pro_684, 
  batch = batch, 
  mod = mod)
# Found10batches
# Adjusting for1covariate(s) or covariate level(s)
# Standardizing Data across genes
# Fitting L/S model and finding priors
# Finding parametric adjustments
# Adjusting the Data

# PCA
pca_noSamplePlate <- prcomp(
  x = t(test_combat_no_SamplePlate), 
  scale. = TRUE)

# Graficar
all(rownames(pca_noSamplePlate$x) == metadata_filtered_684$sampleID)

pdf("pca_no_SamplePlate.pdf")
plot(
  x = pca_noSamplePlate$x[, 1],
  y = pca_noSamplePlate$x[, 2],
  col = as.factor(metadata_filtered_684$batch)
)
dev.off()



# MODELADO DE BATCH COMO BATCH + SENTRIX_ID
batch <- interaction(metadata_filtered_684$batch, metadata_filtered_684$Sentrix_ID)

mod <- model.matrix(~ ceradsc, metadata_filtered_684)

# Remover efecto de lote (batch = SentrixID)
test_combat_no_Batch_noSentrixID <- ComBat(
  dat = m_values_pro_684, 
  batch = batch, 
  mod = mod)
# Using the 'mean only' version of ComBat
# Found128batches
# Note: one batch has only one sample, setting mean.only=TRUE
# Adjusting for1covariate(s) or covariate level(s)
# Error en ComBat(dat = m_values_pro_684, batch = batch, mod = mod): 
#   The covariate is confounded with batch! Remove the covariate and rerun ComBat

# PCA
pca_noBatch_noSentrixID <- prcomp(
  x = t(test_combat_no_Batch_noSentrixID), 
  scale. = TRUE)

# Graficar
all(rownames(pca_noBatch_noSentrixID$x) == metadata_filtered_684$sampleID)

pdf("pca_no_Batch_no_sentrixID.pdf")
plot(
  x = pca_noBatch_noSentrixID$x[, 1],
  y = pca_noBatch_noSentrixID$x[, 2],
  col = as.factor(metadata_filtered_684$batch)
)
dev.off()

###############3#
batch <- interaction(metadata_filtered_684$batch, metadata_filtered_684$Sample_Plate)

mod <- model.matrix(~1, metadata_filtered_684)

# Remover efecto de lote (batch = SentrixID)
test_combat_no_Batch_noSamplePlate_model_Null <- ComBat(
  dat = m_values_pro_684, 
  batch = batch, 
  mod = mod)


#########################
##########################
############################

# Filtrar metadata con diagnóstico Alzheimer (is_AD)
metadata_all_isAD  <- metadata_all %>%
  mutate(is_AD = case_when(
    cogdx == 1 & ceradsc %in% c(3, 4) ~ "control",
    cogdx %in% c(4, 5) & braaksc >= 3 & ceradsc %in% c(1, 2) ~ "AD-NC_SYM",
    TRUE ~ NA_character_
  )) %>% 
  filter(!is.na(is_AD))

# Filtrar valores m con metadata is_AD
sampleIDs <- metadata_all_isAD %>% 
              select(sampleID) %>% 
              unlist() %>% 
              as.vector()

m_values_pro_684_AD <- m_values_pro_684[, colnames(m_values_pro_684) %in% sampleIDs]

dim(m_values_pro_684_AD)
# [1] 335299    357

# Volver a filtrar metadata
sampleIDs <- colnames(m_values_pro_684_AD)

metadata_all_filtered_357 <- metadata_all_isAD %>% 
                              filter(sampleID %in% sampleIDs)

dim(metadata_all_filtered_357)
# [1] 357  50


# Confirmar orden entre metadata y objeto de valores m
all(colnames(m_values_pro_684_AD) == metadata_all_filtered_357$sampleID) #[1] FALSE

# Ordenar metadata en base a objeto valores m
metadata_all_filtered_357 <- metadata_all_filtered_357[match(x = colnames(m_values_pro_684_AD), table = metadata_all_filtered_357$sampleID), ]

# Confirmar orden entre metadata y objeto de valores m
all(colnames(m_values_pro_684_AD) == metadata_all_filtered_357$sampleID) # [1] TRUE

dim(metadata_all_filtered_357)
# [1] 357  50







# library(limma)
# 
# corrected_data <- removeBatchEffect(
#     m_values_pro_684_AD,
#     batch = metadata_all_filtered_357$Sentrix_ID,
#     covariates = model.matrix(~is_AD, metadata_all_filtered_357)[, -1]
#   )
# 
# pca_AD_noSentrixID_corrected_data <- prcomp(
#     x = t(corrected_data),
#     scale. = TRUE
# )
# 
# # Formato para ggplot
# pca_AD_noSentrixID_corrected_data_df <- data.frame(
#   sample = rownames(pca_AD_noSentrixID_corrected_data$x),
#   X = pca_AD_noSentrixID_corrected_data$x[,1],
#   Y = pca_AD_noSentrixID_corrected_data$x[,2]
# )
# 
# pca_pro_var <- pca_AD_noSentrixID_corrected_data$sdev^2
# pca_pro_var_per <- round(pca_pro_var / sum(pca_pro_var) * 100, 1)
# 
# all(pca_AD_noSentrixID_corrected_data_df$sample == metadata_all_filtered_357$sampleID) #[1] TRUE
# 
# # Color batch
# pdf("pca_AD_noSentrixID_corrected_data.pdf")
# pca_AD_noSentrixID_corrected_data_df %>%
#   ggplot(mapping = aes(x = X, y = Y)
#   ) +
#   geom_point() +
#   aes(colour = as.factor(metadata_all_filtered_357$batch)) +
#   scale_color_discrete(name = "batch") +
#   xlab(paste("PC1 - ", pca_pro_var_per[1], "%", sep = "")) +
#   ylab(paste("PC2 - ", pca_pro_var_per[2], "%", sep = "")) +
#   theme_classic() +
#   ggtitle("PCA") +
#   stat_ellipse(geom = "polygon", aes(fill = as.factor(metadata_all_filtered_357$batch)), alpha = 0.2, show.legend = FALSE)
# dev.off()

# PCA SIN CORREGIR EFECTOS DE LOTE

pca_no_corrected_data <- prcomp(
                          x = t(m_values_pro_684_AD),
                     scale. = TRUE
                     )

# Transformar datos para ggplot
pca_no_corrected_data_df <- data.frame(
  sample = rownames(pca_no_corrected_data$x),
  X = pca_no_corrected_data$x[,1],
  Y = pca_no_corrected_data$x[,2]
)

pca_pro_var <- pca_no_corrected_data$sdev^2

pca_pro_var_per <- round(pca_pro_var / sum(pca_pro_var) * 100, 1)

all(pca_no_corrected_data_df$sample == metadata_all_filtered_357$sampleID) # [1] TRUE

# Color batch
pdf("pca_no_corrected_data.pdf")
pca_no_corrected_data_df %>%
  ggplot(mapping = aes(x = X, y = Y)
  ) +
  geom_point() +
  aes(colour = as.factor(metadata_all_filtered_357$Sentrix_ID)) +
  scale_color_discrete(guide = "none") +
  xlab(paste("PC1 - ", pca_pro_var_per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca_pro_var_per[2], "%", sep = "")) +
  theme_classic() +
  ggtitle("PCA (No corregido)") +
  stat_ellipse(geom = "polygon", aes(fill = as.factor(metadata_all_filtered_357$Sentrix_ID)), alpha = 0.2, show.legend = FALSE)
dev.off()


#########################
# CORREGIR EFECTO DE LOTE
#########################
all(metadata_all_filtered_357$sampleID == colnames(m_values_pro_684_AD)) # [1] TRUE

# Corregir efecto de lote conocido
batch <- metadata_all_filtered_357$batch

# Definir modelo de protector
mod <- model.matrix(~as.factor(is_AD), metadata_all_filtered_357)

# Corregir efecto de lote conocido: batch
m_values_pro_357_noBatch <- ComBat(
                              dat = m_values_pro_684_AD,
                            batch = batch, 
                              mod = mod
                            )

# Corregir efecto de lote: sample_plate
batch <- metadata_all_filtered_357$Sample_Plate

# Definir modelo de protector
mod <- model.matrix(~as.factor(is_AD), metadata_all_filtered_357)

m_values_pro_357_noBatch_noSamplePlate <- ComBat(
                              dat = m_values_pro_357_noBatch,
                            batch = batch, 
                              mod = mod
                            )

# Corregir efecto de lote: Desconocido (sva)

# Ajustar modelos
mod <- model.matrix(~is_AD, data = metadata_all_filtered_357)

mod0 <- model.matrix(~1, data = metadata_all_filtered_357)

# SVA para estimar variación oculta
svobj <- sva(
  dat = m_values_pro_357_noBatch_noSamplePlate,
  mod =  mod, 
 mod0 = mod0
) # long time ~ 1.5 hrs

# # Evaluar si SV capturan AD
# apply(svobj$sv, 2, function(sv) summary(lm(sv ~ metadata_all_filtered_357$is_AD))$coefficients[2,4])
# 
# # Evaluar si SV capturan batch
# apply(svobj$sv, 2, function(sv) summary(lm(sv ~ metadata_all_filtered_357$batch))$coefficients[2,4])


# Identificar SV no asociados a AD (se usará para remover en limma)
pvals <- apply(svobj$sv, 2, function(sv)
  summary(lm(sv ~ metadata_all_filtered_357$is_AD))$coefficients[2,4]
)

# AQUI VASSSSS, LUNES 20 ABRIL
svs_noAD <- svobj$sv[, pvals > 0.2]


# Remover variación oculta con limma
data_corrected <- removeBatchEffect(
    m_values_pro_357_noBatch_noSamplePlate,
  covariates = 
)



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
  aes(colour = as.factor(metadata_all_filtered_357$Sentrix_ID)) +
  scale_color_discrete(guide = "none") +
  xlab(paste("PC1 - ", pca_pro_var_per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca_pro_var_per[2], "%", sep = "")) +
  theme_classic() +
  ggtitle("PCA") +
  stat_ellipse(geom = "polygon", aes(fill = as.factor(metadata_all_filtered_357$Sentrix_ID)), alpha = 0.2, show.legend = FALSE)
dev.off()



# BiocManager::install("NOISeq")
# library(NOISeq)
# 
# 
# # read data
# myData <- readData(data = m_values_pro_684_AD, factors = as.data.frame(metadata_all_filtered_357$is_AD))
# 
# 
# myPCA <- dat(input = myData, type = "PCA")
# par(mfrow = c(1,2))
# explo.plot(myPCA, factor = "Tissue")
# explo.plot(myPCA, factor = "batch")




# # batch = metadata_all_filtered_357$batch
# corrected_data <- removeBatchEffect(
#   m_values_pro_684_AD,
#   batch = metadata_all_filtered_357$batch,
#   covariates = model.matrix(~is_AD, metadata_all_filtered_357)[, -1]
# )
# 
# pca_AD_noBatch_corrected_data <- prcomp(
#   x = t(corrected_data),
#   scale. = TRUE)
# 
# # Formato para ggplot
# pca_AD_noBatch_corrected_data_df <- data.frame(
#   sample = rownames(pca_AD_noBatch_corrected_data$x),
#   X = pca_AD_noBatch_corrected_data$x[,1],
#   Y = pca_AD_noBatch_corrected_data$x[,2]
# )
# 
# pca_pro_var <- pca_AD_noBatch_corrected_data$sdev^2
# pca_pro_var_per <- round(pca_pro_var / sum(pca_pro_var) * 100, 1)
# 
# all(pca_AD_noBatch_corrected_data_df$sample == metadata_all_filtered_357$sampleID) #[1] TRUE
# 
# # Color batch
# pdf("pca_AD_noBatch_corrected_data.pdf")
# pca_AD_noBatch_corrected_data_df %>%
#   ggplot(mapping = aes(x = X, y = Y)
#   ) +
#   geom_point() +
#   aes(colour = as.factor(metadata_all_filtered_357$batch)) +
#   scale_color_discrete(name = "Batch") +
#   xlab(paste("PC1 - ", pca_pro_var_per[1], "%", sep = "")) +
#   ylab(paste("PC2 - ", pca_pro_var_per[2], "%", sep = "")) +
#   theme_classic() +
#   ggtitle("PCA") +
#   stat_ellipse(geom = "polygon", aes(fill = as.factor(metadata_all_filtered_357$batch)), alpha = 0.2, show.legend = FALSE)
# dev.off()


# 
# batch <- interaction(metadata_all_filtered_357$batch, metadata_all_filtered_357$Sentrix_ID)
# 
# mod <- NULL
# 
# # Remover efecto de lote
# test_combat_AD_noBatch_noSentrixID<- ComBat(
#   dat = m_values_pro_684_AD,
#   batch = batch,
#   mod = mod)
# 
# # PCA
# pca_AD_noBatch_noSentrixID <- prcomp(
#   x = t(test_combat_AD_noBatch_noSentrixID),
#   scale. = TRUE)
# # Using the 'mean only' version of ComBat
# # Found126batches
# # Note: one batch has only one sample, setting mean.only=TRUE
# # Adjusting for0covariate(s) or covariate level(s)
# # Standardizing Data across genes
# # Error en solve.default(crossprod(design), tcrossprod(t(design), as.matrix(dat))): 
# #   Lapack routine dgesv: system is exactly singular: U[1,1] = 0
#
# NO FUNCIONA


# # Interacción entre batches
# batch <- interaction(metadata_all_filtered_357$batch, metadata_all_filtered_357$Sentrix_ID)
# 
# mod <- model.matrix(~as.factor(is_AD), metadata_all_filtered_357)
# 
# # Remover efecto de lote (batch = SentrixID)
# test_combat_AD_noBatch_noSentrixID<- ComBat(
#   dat = m_values_pro_684_AD, 
#   batch = batch, 
#   mod = mod)
# # Using the 'mean only' version of ComBat
# # Found126batches
# # Note: one batch has only one sample, setting mean.only=TRUE
# # Adjusting for1covariate(s) or covariate level(s)
# # Error en ComBat(dat = m_values_pro_684_AD, batch = batch, mod = mod): 
# #   The covariate is confounded with batch! Remove the covariate and rerun ComBat
#
# NO FUNCIONÓ


# # Interacción entre batches
# batch <- interaction(metadata_all_filtered_357$batch, metadata_all_filtered_357$Sample_Plate)
# 
# mod <- model.matrix(~as.factor(is_AD), metadata_all_filtered_357)
# 
# # Remover efecto de lote (batch = SentrixID)
# test_combat_AD_noBatch_noSamplePlate<- ComBat(
#   dat = m_values_pro_684_AD, 
#   batch = batch, 
#   mod = mod)
# Found20batches
# Adjusting for1covariate(s) or covariate level(s)
# Error en ComBat(dat = m_values_pro_684_AD, batch = batch, mod = mod): 
#   The covariate is confounded with batch! Remove the covariate and rerun ComBat
#
# NO FUNCIONÓ





# # Remover efecto de lote
# batch <- metadata_all_filtered_357$batch
# 
# mod <- model.matrix(~as.factor(is_AD), metadata_all_filtered_357)
# 
# # Remover efecto de lote
# test_combat_AD_noBatch <- ComBat(
#   dat = m_values_pro_684_AD, # En realidad son 357
#   batch = batch, 
#   mod = mod)
# 
# # PCA
# pca_AD_noBatch <- prcomp(
#   x = t(test_combat_AD_noBatch), 
#   scale. = TRUE)
# 
# 
# # Formato para ggplot
# pca_AD_noBatch_df <- data.frame(
#   sample = rownames(pca_AD_noBatch$x),
#   X = pca_AD_noBatch$x[,1],
#   Y = pca_AD_noBatch$x[,2]
# )
# 
# pca_pro_var <- pca_AD_noBatch_df$sdev^2
# pca_pro_var_per <- round(pca_pro_var / sum(pca_pro_var) * 100, 1)
# 
# all(pca_AD_noBatch_df$sample == metadata_all_filtered_357$sampleID) #[1] TRUE
# 
# # Color batch
# pdf("pca_AD_noBatch.pdf")
# pca_AD_noBatch_df %>% 
#   ggplot(mapping = aes(x = X, y = Y)
#   ) +
#   geom_point() +
#   aes(colour = as.factor(metadata_all_filtered_357$batch)) +
#   scale_color_discrete(name = "Batch") +
#   xlab(paste("PC1 - ", pca_pro_var_per[1], "%", sep = "")) +
#   ylab(paste("PC2 - ", pca_pro_var_per[2], "%", sep = "")) +
#   theme_classic() +
#   ggtitle("PCA") +
#   stat_ellipse(geom = "polygon", aes(fill = as.factor(metadata_all_filtered_357$batch)), alpha = 0.2, show.legend = FALSE)
# dev.off()
#
# NO FUNCIONÓ


# # Remover efecto de lote
# batch <- metadata_all_filtered_357$Sentrix_ID
# 
# mod <- model.matrix(~as.factor(is_AD), metadata_all_filtered_357)
# 
# # Remover efecto de lote
# test_combat_AD_noSentrixID <- ComBat(
#   dat = m_values_pro_684_AD, # En realidad son 357
#   batch = batch,
#   mod = mod)
# 
# # PCA
# pca_AD_noSentrixID <- prcomp(
#   x = t(test_combat_AD_noSentrixID),
#   scale. = TRUE)
# 
# # Formato para ggplot
# pca_AD_noSentrixID_df <- data.frame(
#   sample = rownames(pca_AD_noSentrixID$x),
#   X = pca_AD_noSentrixID$x[,1],
#   Y = pca_AD_noSentrixID$x[,2]
# )
# 
# pca_pro_var <- pca_AD_noSentrixID$sdev^2
# pca_pro_var_per <- round(pca_pro_var / sum(pca_pro_var) * 100, 1)
# 
# all(pca_AD_noSentrixID_df$sample == metadata_all_filtered_357$sampleID) #[1] TRUE
# 
# # Color SentrixID
# pdf("pca_AD_noSentrixID.pdf")
# pca_AD_noSentrixID_df %>%
#   ggplot(mapping = aes(x = X, y = Y)
#   ) +
#   geom_point() +
#   aes(colour = as.factor(metadata_all_filtered_357$batch)) +
#   scale_color_discrete(guide = "none") +
#   xlab(paste("PC1 - ", pca_pro_var_per[1], "%", sep = "")) +
#   ylab(paste("PC2 - ", pca_pro_var_per[2], "%", sep = "")) +
#   theme_classic() +
#   ggtitle("PCA") +
#   stat_ellipse(geom = "polygon", aes(fill = as.factor(metadata_all_filtered_357$batch)), alpha = 0.2, show.legend = FALSE)
# dev.off()
# 
# NO FUNCIONÓ




# Con argumento ref.batch en combar // Mismo resultado :(
# #########################################
# # Describir variables al modelo
# batch <- metadata_filtered_684$batch
# 
# mod <- model.matrix(~as.factor(ceradsc), metadata_filtered_684)
# 
# # Remover efecto de lote, modelo ~ cerad
# test_combat_cerad <- ComBat(
#         dat = m_values_pro_684, 
#       batch = batch, 
#         mod = mod, 
#   ref.batch = 1    # Hay 150+ individuos en batch 1
# ) 
# 
# # PCA
# pca_cerad <- prcomp(
#   x = t(test_combat_cerad), 
#   scale. = TRUE)
# 
# # Graficar
# all(rownames(pca_cerad$x) == metadata_filtered_684$sampleID)
# 
# pdf("pca_cerad.pdf")
# plot(
#   x = pca_cerad$x[, 1],
#   y = pca_cerad$x[, 2],
#   col = as.factor(metadata_filtered_684$batch)
# )
# dev.off()
#
# NO FUNCIONÓ


















