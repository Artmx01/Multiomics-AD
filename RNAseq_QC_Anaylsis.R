######################
# QC RNAseq data
######################

# Packages
library(vroom)     # ‘1.7.1’
library(tidyverse) # ‘2.0.0’
library(NOISeq)    # ‘2.56.0’
library(biomaRt)   # ‘2.68.0’
library(EDASeq)    # ‘2.46.0’

# QC RNAseq workflow:
# 1. Filter low counts (CPM > 1)
# 2. Adjust GC bias with EDASeq
# 3. Normalize with TMM method
# 4. Filter genes: only protein-coding, miRNAs, lncRNAs
# 5. Remove batch effect



# 1. Filter low counts (CPM > 1):

# Filter low counts
rnaseq_counts_filtered_215_CPM <- filtered.data(
    dataset = rnaseq_counts_raw_215,
    factor  = factors$is_AD,
       norm = FALSE,
     method = 1,
  cv.cutoff = 100,
        cpm = 1,
      p.adj = "fdr"
)
# Filtering out low count features...
# 15677 features are to be kept for differential expression analysis with filtering method 1

# dim(rnaseq_counts_filtered_215_CPM)
# [1] 15677   215


# 2. Adjust GC bias with EDASeq:

# Build phenoData
phenoData <- factors %>% column_to_rownames(var = "specimenID")

# dim(phenoData)
# [1] 215   3

# Build featureData

# Filter myannot with genes in filtered data (CPM)
myannot_CPM <- myannot %>% filter(ensembl_gene_id %in% rownames(rnaseq_counts_filtered_215_CPM))
myannot_CPM <- myannot_CPM[!duplicated(myannot_CPM$ensembl_gene_id), ]

# dim(myannot_CPM)
# [1] 15624     8

rownames(myannot_CPM) <- NULL
featureData <- myannot_CPM %>% column_to_rownames(var = "ensembl_gene_id")

# dim(featureData)
# [1] 15624     7

# Filtering data counts is needed because some genes are not in featureData
rnaseq_counts_filtered_215_CPM <- rnaseq_counts_filtered_215_CPM[rownames(rnaseq_counts_filtered_215_CPM) %in% rownames(featureData), ]

# dim(rnaseq_counts_filtered_215_CPM)
# [1] 15624   215

# Order featureData based on genes in counts
featureData <- featureData[match(x = rownames(rnaseq_counts_filtered_215_CPM), table = rownames(featureData)), ]

# Check order
identical(rownames(featureData), rownames(rnaseq_counts_filtered_215_CPM))
# [1] TRUE

identical(rownames(phenoData), colnames(rnaseq_counts_filtered_215_CPM))
# [1] TRUE

# Build EDA object
EDA_object <- newSeqExpressionSet(
                     counts = as.matrix(rnaseq_counts_filtered_215_CPM), # matrix input, no data frame
                  phenoData = phenoData,
                featureData = featureData
)

# Adjust GC bias
rnaseq_counts_filtered_215_CPM_noGCbias <- withinLaneNormalization(
              x     = EDA_object,
              y     = "percentage_gene_gc_content",
              which = "full"
)


# 3. Normalize with TMM method

# Create NOISEq object with noGCbias data

# Set names

mylength <- setNames(featureData$length, rownames(featureData))

mybiotype <- setNames(featureData$gene_biotype, rownames(featureData))

mygc <-  setNames(featureData$percentage_gene_gc_content, rownames(featureData))

# Create NOISeq object

noiseqData_filtered_CPM_noGCbias <- NOISeq::readData(
          data    = normCounts(rnaseq_counts_filtered_215_CPM_noGCbias),
          factors = factors,
          length  = mylength,
          biotype = mybiotype,
          gc      = mygc
)

# Normalize data with TMM method

rnaseq_counts_filtered_215_CPM_noGCbias_TMM<- tmm(
            datos = assayData(noiseqData_filtered_CPM_noGCbias)$exprs,
            long = 1000,
            lc = 0,
)

# range(rnaseq_counts_filtered_215_CPM_noGCbias_TMM)
# [1]       0 2016613

# range(rnaseq_counts_raw_215)
# [1]        0 29564642


# 4. Filter genes: only protein-coding, miRNAs, lncRNAs:

# Get genes (features) of interest
keep_protein_coding_miRNAs_lncRNAs <- myannot_CPM %>% filter(gene_biotype %in% c("protein_coding", "miRNA", "lncRNA")) %>% pull(ensembl_gene_id)

# Filter counts with genes of interest (vector)
rnaseq_counts_filtered_215_CPM_noGCbias_TMM_selectedFeatures <- rnaseq_counts_filtered_215_CPM_noGCbias_TMM[rownames(rnaseq_counts_filtered_215_CPM_noGCbias_TMM) %in% keep_protein_coding_miRNAs_lncRNAs, ]

# dim(rnaseq_counts_filtered_215_CPM_noGCbias_TMM_selectedFeatures)
# [1] 15045   215


# 5. Remove batch effect:

# Create NOISeq object with filtered selected features

# Filter myannot: Onlye features of interest
myannot_selectedFeatures <- myannot_CPM %>% filter(gene_biotype %in% c("protein_coding", "miRNA", "lncRNA"))

# dim(myannot_selectedFeatures)
# [1] 15045     8

# Set names

mylength <- setNames(myannot_selectedFeatures$length, myannot_selectedFeatures$ensembl_gene_id)

mybiotype <- setNames(myannot_selectedFeatures$gene_biotype, myannot_selectedFeatures$ensembl_gene_id)

mygc <- setNames(myannot_selectedFeatures$percentage_gene_gc_content, myannot_selectedFeatures$ensembl_gene_id)

# Create NOISeq object

noiseqData_filtered_CPM_noGCbias_TMM_selectedFeatures <- NOISeq::readData(
          data    = rnaseq_counts_filtered_215_CPM_noGCbias_TMM_selectedFeatures,
          factors = factors,
          length  = mylength,
          biotype = mybiotype,
          gc      = mygc
)

# Remove known batch
noiseqData_filtered_215_CPM_noGCbias_TMM_nobatch_selectedFeatures <- ARSyNseq(
    data      = noiseqData_filtered_CPM_noGCbias_TMM_selectedFeatures,
    factor    = "sequencingBatch",
    batch     = TRUE,
    norm      = "n", 
    logtransf = FALSE
)
