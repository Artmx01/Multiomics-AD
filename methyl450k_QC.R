##########################################################################
#                 QC Analysis illumina 450k beadchip
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


