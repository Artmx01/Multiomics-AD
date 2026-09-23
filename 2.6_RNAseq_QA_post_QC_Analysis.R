#####################################
# QA Analysis (postQC): data filtered
#####################################


# Packages
library(vroom)     # ‘1.7.1’
library(tidyverse) # ‘2.0.0’
library(NOISeq)    # ‘2.56.0’
library(biomaRt)   # ‘2.68.0’
library(EDASeq)    # ‘2.46.0’



#      -- Workflow -- 
#
# 1. Check low counts
# 2. Check GC bias
# 3. Check Length bias
# 4. Check RNA compostion bias
# 5. Check meanVar bias
# 6. PCA
#
# Save filtered data



# 1. ------------------------- Check low counts -------------------------

mycountsbio_filtered <- dat(noiseqData_filtered_215_CPM_noGCbias_TMM_selectedFeatures_noBatch,
                   type =  "countsbio",
                   norm = TRUE,
                   factor = NULL)

# Boxplot
pdf("counts_filtered_215.pdf")
explo.plot(mycountsbio_filtered,
           plottype = "boxplot", #type of plot
           samples = 1:50)
dev.off()

# Barplot (Sensitivity Plot)
pdf("counts_filtered_215_barplot.pdf")
explo.plot(mycountsbio_filtered,
           plottype = "barplot", #type of plot
           samples = 1:50)
dev.off()



# 2. ------------------------- Check GC bias -------------------------

 myGCcontent_filtered <- dat(noiseqData_filtered_215_CPM_noGCbias_TMM_selectedFeatures_noBatch,
                   k = 0,            # A feature is considered to be detected if the corresponding number of read counts is > k.
                   type = "GCbias",
                   factor = NULL)

pdf("counts_filtered_215_GCbias.pdf")
explo.plot(myGCcontent_filtered,
           samples = 1:12, # max 12 samples
           toplot = "global")
dev.off()

# GC bias by factor = is_AD

myGCcontent_AD_filtered <- dat(noiseqData_filtered_215_CPM_noGCbias_TMM_selectedFeatures_noBatch,
                   k = 0,
                   type = "GCbias",
                   factor = "is_AD")

pdf("counts_filtered_215_GCbias_AD.pdf")
explo.plot(myGCcontent_AD_filtered,
           samples = NULL,
           toplot = "global")
dev.off()



# 3. ------------------------- Check Length bias -------------------------

 mylengthbias_filtered <- dat(noiseqData_filtered_215_CPM_noGCbias_TMM_selectedFeatures_noBatch,
                    k = 0,
                    type = "lengthbias",
                    factor = NULL)

pdf("counts_filtered_215_lengthBias.pdf")
explo.plot(mylengthbias_filtered,
           samples = 1:12,
           toplot = "global")
dev.off()

# length bias by factor = is_AD

mylengthbias_AD_filtered <- dat(noiseqData_filtered_215_CPM_noGCbias_TMM_selectedFeatures_noBatch,
                    k = 0,
                    type = "lengthbias",
                    factor = "is_AD")

pdf("counts_filtered_215_lengthBias_AD.pdf")
explo.plot(mylengthbias_AD_filtered,
           samples = NULL,
           toplot = "global")
dev.off()



# 4. ------------------------- Check RNA compostion bias -------------------------

rna_comp_bias_filtered <- dat(input = noiseqData_filtered_215_CPM_noGCbias_TMM_selectedFeatures_noBatch, type = "cd", norm = TRUE)

pdf("counts_filtered_215_rna_comp_bias.pdf")
explo.plot(rna_comp_bias_filtered, samples = 1:12)
dev.off()



# 5. ------------------------- Check meanVar bias -------------------------

# Average gene expression
avg_gene_exps_filtered_215 <- rowMeans(exprs(noiseqData_filtered_215_CPM_noGCbias_TMM_selectedFeatures_noBatch))

# Variance per gene
gene_vars_filtered_215 <- apply(exprs(noiseqData_filtered_215_CPM_noGCbias_TMM_selectedFeatures_noBatch), 1, var)

# Plot mean-variance relation (log-log scale)
pdf("meanVar_filtered_215.pdf")
plot(x = avg_gene_exps_filtered_215, y = gene_vars_filtered_215, log = "xy", main = "Mean-Variance relation", xlab = "Mean expression (log scale)", ylab = "Variance (log scale)")
dev.off()

# Scaling is needed (performed in sgcca function (wrapper.sgcca))



# 6. ------------------------- PCA -------------------------

# Get PCA
pca_rnaseq_filtered_215 <- prcomp(
  x      = t(exprs(noiseqData_filtered_215_CPM_noGCbias_TMM_selectedFeatures_noBatch)), 
  center = TRUE,
  scale. = FALSE
)

# PCA dataframe (for plotting)
pca_rnaseq_filtered_215_df <- data.frame(
  sample = rownames(pca_rnaseq_filtered_215$x),
  X      = pca_rnaseq_filtered_215$x[,1],
  Y      = pca_rnaseq_filtered_215$x[,2]
)

# Get % variance per PC
pca_var_rnaseq_filtered_215 <- pca_rnaseq_filtered_215$sdev^2

pca_var_rnaseq_filtered_215_per <- round(pca_var_rnaseq_filtered_215 / sum(pca_var_rnaseq_filtered_215) * 100 , 1)

identical(pca_rnaseq_filtered_215_df$sample, RNA_seq_metadata_filteredQC_DLPFC_215$specimenID)
# [1] TRUE

# Plot PCA by ADcontrol
pdf("pca_rnaseq_filtered_215_isAD.pdf")
  pca_rnaseq_filtered_215_df %>%
 ggplot(mapping = aes(x = X, y = Y)
      ) +
      geom_point() +
      aes(colour = as.factor(RNA_seq_metadata_filteredQC_DLPFC_215$is_AD)) +
      scale_color_discrete(name = "AD/Control") +
      xlab(paste("PC1 - ", pca_var_rnaseq_filtered_215_per[1], "%", sep = "")) +
      ylab(paste("PC2 - ", pca_var_rnaseq_filtered_215_per[2], "%", sep = "")) +
      theme_classic() +
      ggtitle("PCA") +
      stat_ellipse(geom = "polygon", aes(fill = as.factor(RNA_seq_metadata_filteredQC_DLPFC_215$is_AD)), alpha = 0.2, show.legend = FALSE)
dev.off()

# PCA (Colour = sequencingBatch)
pdf("pca_rnaseq_filtered_215_batch.pdf")
  pca_rnaseq_filtered_215_df %>%
 ggplot(mapping = aes(x = X, y = Y)
      ) +
      geom_point() +
      aes(colour = as.factor(RNA_seq_metadata_filteredQC_DLPFC_215$sequencingBatch)) +
      scale_color_discrete(name = "Batch") +
      xlab(paste("PC1 - ", pca_var_rnaseq_filtered_215_per[1], "%", sep = "")) +
      ylab(paste("PC2 - ", pca_var_rnaseq_filtered_215_per[2], "%", sep = "")) +
      theme_classic() +
      ggtitle("PCA") +
      stat_ellipse(geom = "polygon", aes(fill = as.factor(RNA_seq_metadata_filteredQC_DLPFC_215$sequencingBatch)), alpha = 0.2, show.legend = FALSE)
dev.off()


# ------------------------- QA Analysis finished (post QC) -------------------------

# Save filtered data
saveRDS(object = exprs(noiseqData_filtered_215_CPM_noGCbias_TMM_selectedFeatures_noBatch), file = "ROSMAP_RNAseq_counts_filtered_215.rds")

# Save metadata filtered
vroom_write(RNA_seq_metadata_filteredQC_DLPFC_215 , file = "ROSMAP_RNAseq_metadata_filtered_215.txt")

# Save annotation features (of genes)
vroom_write(myannot_selectedFeatures, file = "ROSMAP_RNAseq_annotationFeatures_filtered_215.txt")
