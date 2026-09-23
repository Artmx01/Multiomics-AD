######################################
# RNAseq QA ANALYSIS (with raw counts)
######################################


# Packages
library(vroom)     # ‘1.7.1’
library(tidyverse) # ‘2.0.0’
library(NOISeq)    # ‘2.56.0’
library(biomaRt)   # ‘2.68.0’
library(EDASeq)    # ‘2.46.0’



#                  -- Worflow --
#
# 1.----- Metadata -----
# 1.1 Load metadata
#
# 2.----- Data -----
# 2.1 Load data
# 2.2 Remove duplicated rows (median expression)
# 2.3 Filter data (3 omics & AD/control)
#
# 3.----- QA Analysis -----
# 3.1 Mart Annotation 
# 3.2 Create NOISeq object
# 3.3 Check low counts
# 3.4 GC bias
# 3.5 Length bias
# 3.6 RNA composition bias
# 3.7 meanVar plot
# 3.8 PCA
#
# QA Analysis finished



# 1. ------------------------- Metadata -------------------------

# Load metadata:
RNA_seq_metadata_filteredQC_DLPFC <- vroom(file = "/STORAGE/csbig/multiomics-Arturo/mRNA/RNA_seq_metadata_filteredQC_DLPFC.txt")

# Filter metadata (subjects with 3 omics & AD/Control):
RNA_seq_metadata_filteredQC_DLPFC_215 <- RNA_seq_metadata_filteredQC_DLPFC %>%
    filter(individualID %in% metadata_filtered_isAD$individualID) %>% 
    distinct(individualID, .keep_all = TRUE) %>% 
    mutate(is_AD = case_when(
      cogdx == 1 & ceradsc %in% c(3, 4) ~ "control",
      cogdx %in% c(4, 5) & braaksc >= 3 & ceradsc %in% c(1, 2) ~ "AD-NC_SYM",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(is_AD))
    
# dim(RNA_seq_metadata_filteredQC_DLPFC_215)
# [1] 215  42



# 2. ------------------------- Data -------------------------

# 2.1 Load data:

rnaseq_counts_raw <- readRDS(file = "/STORAGE/csbig/multiomics-Arturo/mRNA/ROSMAP_RNAseq_rawcounts_DLPFC.rds")


# 2.2 Remove duplicated rows (with his median expression):

# As in: 
# https://github.com/paulinapglz99/Multinetwork-topologic-analysis/blob/main/00_preprocessing/2.pre-pro-mRNA_ARSyn.R


# strip transcript or gene version numbers from identifiers 
rnaseq_counts_raw <- rnaseq_counts_raw %>%
    mutate(feature = str_remove(feature, "\\..*$"))

# Get repeated features
repeated_features <- rnaseq_counts_raw %>%
    group_by(feature) %>%
    filter(n() > 1) %>%
    distinct(feature) %>%
    pull(feature)

# length(repeated_features)
# [1] 45

# View duplicated rows
repeated_rows <- rnaseq_counts_raw[rnaseq_counts_raw$feature %in% repeated_features, ]

# Sort and calculate the median of duplicate values
repeated_rows <- repeated_rows[order(repeated_rows$feature),]

repeated_rows <- repeated_rows %>%
    group_by(feature) %>%
    summarize(across(everything(), \(x) median(x, na.rm = TRUE)))

# Remove duplicated rows and add rows with the calculated median
rnaseq_counts_raw <- rnaseq_counts_raw %>% filter(!feature %in% repeated_rows$feature)

rnaseq_counts_raw <- bind_rows(rnaseq_counts_raw, repeated_rows)

# Set feature column as rownames
rnaseq_counts_raw <- column_to_rownames(rnaseq_counts_raw, var = "feature")

# dim(rnaseq_counts_raw)
# [1] 60562  1141


# 2.3 Filter data (3 omics & AD/control)

# Get subjects
subjects <- RNA_seq_metadata_filteredQC_DLPFC_215$specimenID

# Filter raw data (subjects with 3 omics & AD/Control)
rnaseq_counts_raw_215 <- rnaseq_counts_raw[,subjects]

# dim(rnaseq_counts_raw_215)
# [1] 60562   215



# 3. ------------------------- QA Analysis -------------------------

# 3.1  Mart Annotation:

# Get mart
mart <- useEnsembl("ensembl", dataset="hsapiens_gene_ensembl", version = 110)

# Get annotation features
myannot <- getBM(attributes = c("ensembl_gene_id", "chromosome_name",
                                "percentage_gene_gc_content", "gene_biotype",
                                "start_position","end_position","hgnc_symbol"),
                 filters = "ensembl_gene_id",
                 values  =  rownames(rnaseq_counts_raw_215),
                 mart    = mart)


# 3.2 Create NOISeq object:

# Set factors object (necessary in NOISeq object)
factors <- RNA_seq_metadata_filteredQC_DLPFC_215 %>% 
            dplyr::select(specimenID, is_AD, sequencingBatch, libraryBatch)

factors <- as.data.frame(factors)

# Check order btw data & factors
identical(colnames(rnaseq_counts_raw_215), factors$specimenID)
# [1] TRUE

# Set attributes (length, GC, biotype, chr)

myannot$length <- abs(myannot$end_position - myannot$start_position + 1)

mylength <- setNames(myannot$length, myannot$ensembl_gene_id)

mygc <- setNames(myannot$percentage_gene_gc_content, myannot$ensembl_gene_id)

mybiotype <-setNames(myannot$gene_biotype, myannot$ensembl_gene_id)

## Create NOISeq object with attributes
noiseqData_raw <- NOISeq::readData(data = rnaseq_counts_raw_215,
                               factors = factors,
                               gc = mygc,
                               biotype = mybiotype,
                               length =  mylength)


# Diagnostic Plots:

# 3.3 Check low counts:

mycountsbio_raw <- dat(noiseqData_raw,
                   type =  "countsbio",
                   norm = F,
                   factor = NULL)

# boxplot
pdf("Counts_raw_215.pdf")
explo.plot(mycountsbio_raw,
           plottype = "boxplot", #type of plot
           samples = 1:50)
dev.off()

# Sensitivity plot (barplot)
pdf("counts_raw_215_barplot.pdf")
explo.plot(mycountsbio_raw,
           plottype = "barplot",
           samples = 1:50)
dev.off()


# 3.4 GC bias:

myGCcontent_raw <- dat(noiseqData_raw,
                   k = 0,            # A feature is considered to be detected if the corresponding number of read counts is > k.
                   type = "GCbias",
                   factor = NULL)

pdf("counts_raw_215_GCbias.pdf")
explo.plot(myGCcontent_raw,
           samples = 1:12, # max 12 samples
           toplot = "global")
dev.off()

# GC bias by AD/Control
myGCcontent_AD_raw <- dat(noiseqData_raw,
                   k = 0,
                   type = "GCbias",
                   factor = "is_AD")

pdf("counts_raw_215_GCbias_AD.pdf")
explo.plot(myGCcontent_AD_raw,
           samples = NULL,
           toplot = "global")
dev.off()


# 3.5 Length bias:

mylengthbias_raw <- dat(noiseqData_raw,
                    k = 0,
                    type = "lengthbias",
                    factor = NULL)

pdf("count_raw_215_lengthBias.pdf")
explo.plot(mylengthbias_raw,
           samples = 1:12,
           toplot = "global")
dev.off()

# Length bias by AD/Control
mylengthbias_AD_raw <- dat(noiseqData_raw,
                    k = 0,
                    type = "lengthbias",
                    factor = "is_AD")

pdf("count_raw_215_lengthBias_AD.pdf")
explo.plot(mylengthbias_AD_raw,
           samples = NULL,
           toplot = "global")
dev.off()


# 3.6 RNA composition bias:

rna_comp_bias <- dat(input = noiseqData_raw, type = "cd", norm = FALSE)

# Plot RNA composition bias
pdf("counts_raw_215_rna_comp_bias.pdf")
explo.plot(rna_comp_bias, samples = 1:12)
dev.off()


# 3.7 meanVar plot:

# Average gene expresión
avg_gene_exps_raw_215 <- rowMeans(rnaseq_counts_raw_215)

# Variance per gene
gene_vars_raw_215 <- apply(rnaseq_counts_raw_215, 1, var)

# Plot mean variance relation (log-log scale)
pdf("meanVar_raw_215.pdf")
plot(
    x    = avg_gene_exps_raw_215, 
    y    = gene_vars_raw_215, 
    log  = "xy", 
    main = "Mean-Variance relation", 
    xlab = "Mean expression (log scale)", 
    ylab = "Variance (log scale)")
dev.off()


# 3.8 PCA:

# get PCA
pca_rnaseq_raw_215 <- prcomp(x = t(rnaseq_counts_raw_215), center = TRUE, scale. = FALSE)

# PCA (data frame, for plotting)
pca_rnaseq_raw_215_df <- data.frame(
  sample = rownames(pca_rnaseq_raw_215$x),
  X = pca_rnaseq_raw_215$x[,1],
  Y = pca_rnaseq_raw_215$x[,2]
)

# Get variance for each PCs
pca_var_rnaseq_raw_215 <- pca_rnaseq_raw_215$sdev^2

# Get % variance for each PCs
pca_var_rnaseq_raw_215_per <- round(pca_var_rnaseq_raw_215 / sum(pca_var_rnaseq_raw_215) * 100 , 1)

# Check order between samples and metadata (important for colouring)
identical(pca_rnaseq_raw_215_df$sample, RNA_seq_metadata_filteredQC_DLPFC_215$specimenID)
# [1] TRUE

# Plot
pdf("pca_rnaseq_raw_215_isAD.pdf")
  pca_rnaseq_raw_215_df %>%
 ggplot(mapping = aes(x = X, y = Y)
      ) +
      geom_point() +
      aes(colour = as.factor(RNA_seq_metadata_filteredQC_DLPFC_215$is_AD)) +
      scale_color_discrete(name = "AD/Control") +
      xlab(paste("PC1 - ", pca_var_rnaseq_raw_215_per[1], "%", sep = "")) +
      ylab(paste("PC2 - ", pca_var_rnaseq_raw_215_per[2], "%", sep = "")) +
      theme_classic() +
      ggtitle("PCA") +
      stat_ellipse(geom = "polygon", aes(fill = as.factor(RNA_seq_metadata_filteredQC_DLPFC_215$is_AD)), alpha = 0.2, show.legend = FALSE)
dev.off()

# PCA (color = sequencingBatch)

# ------------------------ QA Analysis finished -------------------------
