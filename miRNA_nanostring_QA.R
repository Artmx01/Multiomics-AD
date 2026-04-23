##########################################################################
#                 QA Analysis miRNA nanostring nCounter
##########################################################################
# 13 abril 2026
# Arturo, BM

# BiocManager::install("NanoStringNCTools")
# install.packages('NanoStringNorm')
# devtools::install_github("uclahs-cds/package-NanoStringNorm")
# # Cargar librerías
# library(NanoStringNorm)
# library(NanoStringNCTools)


# Contruir metadata_all_miRNA







# Lectura de datos miRNA (gct)
# gct format: matrix con rows (miRNAs) x columns (samples). Contiene las abundancias absolutas de miRNA
# por muestra.


# LECTURA DE ARCHIVOS RAW (RCC FILES; REPORTER CODE COUNT FILES)
#Lectura de archivos
miRNA_data <- read.delim(
        file = "/STORAGE/csbig/multiomics-Arturo/miRNA_data/ROSMAP_arraymiRNA.gct",
        skip = 2,
      header = TRUE
)

# Dimensiones del data frame
dim(miRNA_data)
# [1] 309 704

# Asignar nombre a las filas (miRNAS) utilizando la columna Genes
rownames(miRNA_data) <- miRNA_data$Genes

#Eliminar las primeras 2 columnas (solo mantener samples)
miRNA_data <- miRNA_data[,-c(1,2)]

#Confirmar filtrado
dim(miRNA_data) 
# [1] 309 702 // 309 miRNAs x 702 individuos

# Parece ser que ya están preprocesados
# https://www.synapse.org/Synapse:syn3387325





# Filtrar metadata
metadata_all_miRNA_isAD <- metadata_all_miRNA %>% 
                            mutate(is_AD = case_when(
                              cogdx == 1 & ceradsc %in% c(3, 4) ~ "control",
                              cogdx %in% c(4, 5) & braaksc >= 3 & ceradsc %in% c(1, 2) ~ "AD-NC_SYM",
                              TRUE ~ NA_character_
                            )) %>% 
                              filter(!is.na(is_AD)) %>% 
                              filter(!is.na(mirna_id...2))

dim(metadata_all_miRNA_isAD)
# [1] 267  51


# miRNA_samples <- colnames(miRNA_data)
# 
# miRNA_metadata_filtered_525 <- miRNA_metadata %>% 
#   filter(mirna_id...2 %in% miRNA_samples)
# 
# dim(miRNA_metadata_filtered_525)
# # [1] 525  14
# 
# # Agregar columna is_AD a metadata
# View(metadata_all_filtered_357)
# View(miRNA_metadata_filtered_525)
# # Usaremos metadata_all_filtered_357 (450k beadchip)



# 
# # Filtrar miRNA data, solo sujetos en metadata_filtered
# keep <- miRNA_metadata_filtered_525 %>% 
#           select(mirna_id...2) %>% 
#           unlist() %>% 
#           as.vector()
# 
# # ¿Cuántos individuos son?
# length(keep)
# # [1] 525
# 
# # Filtrar miRNA data
# miRNA_data_filtered_525 <- miRNA_data[,keep]
# 
# # Confirmar filtrado
# dim(miRNA_data_filtered_525)
# # [1] 309 525


# # Filtrar miRNA data, solo sujetos en metadata_all_miRNA_isAD

keep <- metadata_all_miRNA_isAD %>% 
  select(mirna_id...2) %>% 
  unlist() %>% 
  as.vector()

length(keep) # [1] 267 // Algunos no contienen info en la columna mirna_id...2

####################################
#IMPORTANTE REALIZAAR:
# Lo más adecuado sería modficar colnames para cambiar mirna_id...2 a individualID/biospecimen!!!!!
# Así perderíamos menos muestras
####################################

# Filtrar miRNA data
miRNA_data_filtered_266 <- miRNA_data[,colnames(miRNA_data) %in% keep]

# Confirmar filtrado
dim(miRNA_data_filtered_266)
# [1] 309 266


# QA ANALYSIS
############################

# ¿Existen NAs?
sum(is.na(miRNA_data_filtered_266))
# [1] 0 

# ¿Cuál es el rango de los datos?
range(miRNA_data_filtered_266)
# [1]  2.444961 16.501511 // Indica preprocesamiento previo (log2)

# Resumen de los datos
summary(miRNA_data_filtered_266)

# Boxplot de los datos (50 muestras)
pdf("boxplot_miRNA.pdf")
boxplot(x = miRNA_data_filtered_266[,1:85], outline = FALSE)
dev.off()


# PCA

miRNA_pca_266 <- prcomp(
                x = t(miRNA_data_filtered_266), 
           scale. = TRUE
           )

# Tranformar datos para ggplot
miRNA_pca_266_df <- data.frame(
           sample = rownames(miRNA_pca_266$x),
           X = miRNA_pca_266$x[,1],
           Y = miRNA_pca_266$x[,2]
           )

# Extraer variabilidad de cada PC
pca_var <- miRNA_pca_266$sdev^2

# Variabilidad en %
pca_var_per <- round(pca_var / sum(pca_var) * 100, 1)


#Volver a filtrar metadata
samples <- rownames(miRNA_pca_266$x)


metadata_all_miRNA_isAD_266 <- metadata_all_miRNA_isAD %>%
  filter(mirna_id...2 %in% samples)

# Confirmar filtrado
dim(metadata_all_miRNA_isAD_266)

# Evaluar orden
all(miRNA_pca_266_df$sample == metadata_all_miRNA_isAD_266$mirna_id...2) # [1] FALSE



pdf("miRNA_pca.pdf")
miRNA_pca_525_df %>%
  ggplot(mapping = aes(x = X, y = Y)
  ) +
  geom_point() +
  aes(colour = as.factor(miRNA_metadata_filtered_525$plate)) +
  scale_color_discrete(name = "Plate") +
  xlab(paste("PC1 - ", pca_var_per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca_var_per[2], "%", sep = "")) +
  theme_classic() +
  ggtitle("PCA") +
  stat_ellipse(geom = "polygon", aes(fill = as.factor(miRNA_metadata_filtered_525$plate)), alpha = 0.2, show.legend = FALSE)
dev.off()

# silhouette Score (plate)
miRNA.silhouette <- batch_sil(miRNA_pca_525, miRNA_metadata_filtered_525$plate)
miRNA.silhouette
# [1] -0.06742445













