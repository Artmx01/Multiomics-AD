####################################
# Fit SGCCA for AD data & stability
####################################


# libraries
library(RGCCA)        # ‘3.0.3’
library(RhpcBLASctl)  # ‘0.23.42’
library(tidyverse)    # ‘2.0.0’









# Fit SGCCA model for AD dataset:
fit_sgcca_AD <- rgcca(
  blocks      = sgcca_AD_data,
  connection  = design,
  ncomp       = ncomp_permut_50$best_params,
  sparsity    = penalties_permut_500$best_params,
  scheme      = "centroid",
  scale       = TRUE,
  scale_block = "lambda1",
  method      = "sgcca"
)

print(fit_sgcca_AD)
# Fitted SGCCA model.
# The algorithm converged to a stationnary point after 41 iterations.

summary(fit_sgcca_AD)
# Call: method='sgcca', superblock=FALSE, scale=TRUE, scale_block='lambda1',
# init='svd', bias=TRUE, tol=1e-08, NA_method='na.ignore', ncomp=c(1,1,1),
# response=NULL, comp_orth=TRUE
# There are J = 3 blocks.
# The design matrix is:
#       methyl rnaseq mirna
# methyl      0      1     1
# rnaseq      1      0     1
# mirna       1      1     0

# The centroid scheme is used.
# Sum_{j,k} c_jk g(cov(X_j a_j, X_k a_k) = 0.1917

# The sparsity parameter used for methyl is: 0.08490206 (with 4453
# variables selected)
# The sparsity parameter used for rnaseq is: 0.09080668 (with 208
# variables selected)
# The sparsity parameter used for mirna is: 0.1354807 (with 8 variables selected)


# Model explains > 50% variance between blocks (AVE_inner):
fit_sgcca_AD$AVE
# $AVE_X
# $AVE_X$methyl
# [1] 0.01160635

# $AVE_X$rnaseq
# [1] 0.2045678

# $AVE_X$mirna
# [1] 0.03412276


# $AVE_X_cor
# $AVE_X_cor$methyl
# [1] 0.01160635

# $AVE_X_cor$rnaseq
# [1] 0.2045678

# $AVE_X_cor$mirna
# [1] 0.03412276


# $AVE_outer
# [1] 0.02500157

# $AVE_inner
# [1] 0.5568738



# Select most stable variables

boot_sgcca_AD <- rgcca_bootstrap(
                  fit_sgcca_AD,
                  n_boot             = 1000,
                  n_cores            = 50,
                  balanced           = TRUE,
                  keep_all_variables = FALSE,
                  verbose            = TRUE
)




features_freq_AD_df <-  boot_sgcca_AD$bootstrap %>% 
                            filter(type == "weights") %>%
                            group_by(block, var) %>%
                            summarise( 
                               freq   = mean(value != 0),
                              .groups = "drop"
                            )
                            
stable_features_AD_df <- boot_sgcca_AD$stats %>%
                           filter(type == "weights") %>%
                           left_join(features_freq_AD_df, by = c("block", "var")) %>%
                           group_by(block) %>%
                           mutate (pval.adjust = p.adjust(pval, method = "BH")) %>%
                           ungroup() %>%
                           filter(freq >= 0.8, pval.adjust < 0.05)

dim(stable_features_AD_df)
# [1] 4663   16 // 4663 features (CpGs, mRNAs, miRNAs) from which:

table(stable_features_AD_df$block)
# methyl  mirna rnaseq
#  4447      8    208











  

# EN OTRO SCRIPT:


# Fit SGCCA model for control dataset:
fit_sgcca_control <- rgcca(
  blocks      = sgcca_control_data,
  connection  = design,
  ncomp       = ncomp_permut_50$best_params,
  sparsity    = penalties_permut_500$best_params,
  scheme      = "centroid",
  scale       = TRUE,
  scale_block = "lambda1",
  method      = "sgcca"
)

print(fit_sgcca_control)
# Fitted SGCCA model.
# The algorithm converged to a stationnary point after 31 iterations.

summary(fit_sgcca_control)
# Call: method='sgcca', superblock=FALSE, scale=TRUE, scale_block='lambda1',
# init='svd', bias=TRUE, tol=1e-08, NA_method='na.ignore', ncomp=c(1,1,1),
# response=NULL, comp_orth=TRUE
# There are J = 3 blocks.
# The design matrix is:
#       methyl rnaseq mirna
# methyl      0      1     1
# rnaseq      1      0     1
# mirna       1      1     0

# The centroid scheme is used.
# Sum_{j,k} c_jk g(cov(X_j a_j, X_k a_k) = 0.1983

# The sparsity parameter used for methyl is: 0.08490206 (with 4523
# variables selected)
# The sparsity parameter used for rnaseq is: 0.09080668 (with 206
# variables selected)
# The sparsity parameter used for mirna is: 0.1354807 (with 10
# variables selected)


# Model explains > 50% variance between blocks (AVE_inner):
fit_sgcca_control$AVE
# $AVE_X
# $AVE_X$methyl
# [1] 0.01670228

# $AVE_X$rnaseq
# [1] 0.2274933

# $AVE_X$mirna
# [1] 0.1201115


# $AVE_X_cor
# $AVE_X_cor$methyl
# [1] 0.01670228

# $AVE_X_cor$rnaseq
# [1] 0.2274933

# $AVE_X_cor$mirna
# [1] 0.1201115


# $AVE_outer
# [1] 0.04366314

# $AVE_inner
# [1] 0.5341872


boot_sgcca_control <- rgcca_bootstrap(
                  fit_sgcca_control,
                  n_boot             = 1000,
                  n_cores            = 50,
                  balanced           = TRUE,
                  keep_all_variables = FALSE,
                  verbose            = TRUE
)

features_freq_control_df <-  boot_sgcca_control$bootstrap %>% 
                            filter(type == "weights") %>%
                            group_by(block, var) %>%
                            summarise( 
                               freq   = mean(value != 0),
                              .groups = "drop"
                            )

stable_features_control_df <- boot_sgcca_control$stats %>%
                           filter(type == "weights") %>%
                           left_join(features_freq_control_df, by = c("block", "var")) %>%
                           group_by(block) %>%
                           mutate (pval.adjust = p.adjust(pval, method = "BH")) %>%
                           ungroup() %>%
                           filter(freq >= 0.8, pval.adjust < 0.05)


dim(stable_features_control_df)
# [1] 4508   16 // 4508 features, from which:

table(stable_features_control_df$block)
# methyl  mirna rnaseq
#  4292     10    206



# Get stable features data (AD) for network reconstruction

CpGs_stable_AD <- stable_features_AD_df %>% filter(block == "methyl") %>% pull(var)
rnaseq_stable_AD <- stable_features_AD_df %>% filter(block == "rnaseq") %>% pull(var)
mirna_stable_AD <- stable_features_AD_df %>% filter(block == "mirna") %>% pull(var)

stable_features_AD_data <- list(
                            methyl = methyl_AD_data[, CpGs_stable_AD],
                            rnaseq = rnaseq_AD_data[, rnaseq_stable_AD],
                            mirna  = mirna_AD_data[, mirna_stable_AD]
                           )



# Fit sgcca with stable features (AD)
fit_sgcca_AD_stab_features <- rgcca(
  blocks      = stable_features_AD_data,
  connection  = design,
  ncomp       = ncomp_permut_50$best_params,
  sparsity    = 1, # Keep all variables
  scheme      = "centroid",
  scale       = TRUE,
  scale_block = "lambda1",
  method      = "sgcca"
)

# Extract data blocks from fitted model:
# These blocks are preprocessed according to the values of scale/scale_block/NA_method
stable_features_AD_data_Networks <- fit_sgcca_AD_stab_features$blocks

# transpose blocks
stable_features_AD_data_Networks <- lapply(stable_features_AD_data_Networks, t)

# 1. Get methyl block and save it as a tsv format (with rownames as a column)
methyl_AD_data_networks <- stable_features_AD_data_Networks$methyl

methyl_AD_data_networks <- data.frame(
                                 gene = rownames(methyl_AD_data_networks),
                                        as.data.frame(methyl_AD_data_networks),
                            row.names = NULL
                          )
# Quality Check
identical(rownames(stable_features_AD_data_Networks$methyl), methyl_AD_data_networks$gene)
# [1] TRUE

# Save methyl AD data
write_tsv(x = methyl_AD_data_networks, file = "methyl_AD_data_networks.tsv")


# 2. Get rnaseq block and save it as a tsv format (with rownames as a column)

rnaseq_AD_data_networks <- stable_features_AD_data_Networks$rnaseq

rnaseq_AD_data_networks <- data.frame(
                                 gene = rownames(rnaseq_AD_data_networks),
                                        as.data.frame(rnaseq_AD_data_networks),
                            row.names = NULL
                          )

identical(rownames(stable_features_AD_data_Networks$rnaseq), rnaseq_AD_data_networks$gene)
# [1] TRUE

# Save rnaseq AD data
write_tsv(x = rnaseq_AD_data_networks, file = "rnaseq_AD_data_networks.tsv")


# 3. Get mirna block and save it as a tsv format (with rownames as a column)

mirna_AD_data_networks <- stable_features_AD_data_Networks$mirna

mirna_AD_data_networks <- data.frame(
                                 gene = rownames(mirna_AD_data_networks),
                                        as.data.frame(mirna_AD_data_networks),
                            row.names = NULL
                          )

identical(rownames(stable_features_AD_data_Networks$mirna), mirna_AD_data_networks$gene)
# [1] TRUE

# Save mirna AD data
write_tsv(x = mirna_AD_data_networks, file = "mirna_AD_data_networks.tsv")

# Save all AD data (stable features)
saveRDS(object = stable_features_AD_data_Networks, file = "stable_features_AD_data_Networks.rds")




# Get stable features data (control) for network reconstruction

CpGs_stable_control <- stable_features_control_df %>% filter(block == "methyl") %>% pull(var)
rnaseq_stable_control <- stable_features_control_df %>% filter(block == "rnaseq") %>% pull(var)
mirna_stable_control <- stable_features_control_df %>% filter(block == "mirna") %>% pull(var)

stable_features_control_data <- list(
                                  methyl = methyl_control_data[, CpGs_stable_control],
                                  rnaseq = rnaseq_control_data[, rnaseq_stable_control],
                                  mirna  = mirna_control_data[, mirna_stable_control]
                                )


# Fit sgcca with stable features (control)
fit_sgcca_control_stab_features <- rgcca(
  blocks      = stable_features_control_data,
  connection  = design,
  ncomp       = ncomp_permut_50$best_params,
  sparsity    = 1, # Keep all variables
  scheme      = "centroid",
  scale       = TRUE,
  scale_block = "lambda1",
  method      = "sgcca"
)

# Extract data blocks from fitted model:
# These blocks are preprocessed according to the values of scale/scale_block/NA_method
stable_features_control_data_Networks <- fit_sgcca_control_stab_features$blocks

# transpose blocks
stable_features_control_data_Networks <- lapply(stable_features_control_data_Networks, t)

# 1. Get methyl block and save it as a tsv format (with rownames as a column)
methyl_control_data_networks <- stable_features_control_data_Networks$methyl

methyl_control_data_networks <- data.frame(
                                 gene = rownames(methyl_control_data_networks),
                                        as.data.frame(methyl_control_data_networks),
                            row.names = NULL
                          )

identical(rownames(stable_features_control_data_Networks$methyl), methyl_control_data_networks$gene)
# [1] TRUE

# Save methyl control data
write_tsv(x = methyl_control_data_networks, file = "methyl_control_data_networks.tsv")


# 2. Get rnaseq block and save it as a tsv format (with rownames as a column)
rnaseq_control_data_networks <- stable_features_control_data_Networks$rnaseq

rnaseq_control_data_networks <- data.frame(
                                 gene = rownames(rnaseq_control_data_networks),
                                        as.data.frame(rnaseq_control_data_networks),
                            row.names = NULL
                          )

identical(rownames(stable_features_control_data_Networks$rnaseq), rnaseq_control_data_networks$gene)
# [1] TRUE

# Save rnaseq control data
write_tsv(x = rnaseq_control_data_networks, file = "rnaseq_control_data_networks.tsv")


# 3. Get mirna block and save it as a tsv format (with rownames as a column)
mirna_control_data_networks <- stable_features_control_data_Networks$mirna

mirna_control_data_networks <- data.frame(
                                 gene = rownames(mirna_control_data_networks),
                                        as.data.frame(mirna_control_data_networks),
                            row.names = NULL
                          )

identical(rownames(stable_features_control_data_Networks$mirna), mirna_control_data_networks$gene)
# [1] TRUE

# Save rnaseq control data
write_tsv(x = mirna_control_data_networks, file = "mirna_control_data_networks.tsv")


# Save all control data (stable features)
saveRDS(object = stable_features_control_data_Networks, file = "stable_features_control_data_Networks.rds")






# Venn Diagram of features selected by SGCCA in AD & control:

library(ggVennDiagram) # ‘1.5.7’
library(ggplot2)       # ‘4.0.3’


# Features stable in AD
stable_features_AD_names <- c(CpGs_stable_AD, rnaseq_stable_AD, mirna_stable_AD)

# Features stable in control
stable_features_control_names <- c(CpGs_stable_control, rnaseq_stable_control, mirna_stable_control)

# Venn diagram data
Venn_data <- list(AD = stable_features_AD_names, noAD = stable_features_control_names)

# Plot and save Venn diagram
pdf("VennDiagram_AD_noAD_features_selected.pdf")
ggVennDiagram(Venn_data)
dev.off()





###############################
# Functional Enrichment script
###############################

library(clusterProfiler) # ‘4.20.3’
library(biomaRt)         # ‘2.68.0’
library(vroom)           # ‘1.7.1’
library(org.Hs.eg.db)    # ‘3.23.1’
library(tidyverse)       # ‘2.0.0’






# Load  AD data (stable features from SGCCA):
stable_features_AD_data <- readRDS(file = "../networks/stable_features_AD_data_Networks.rds")

# Load control data (stable features from SGCCA):
stable_features_control_data <- readRDS(file = "../networks/stable_features_control_data_Networks.rds")

#Load methyl array metadata
array_meth_metadata <- vroom(file = "/STORAGE/csbig/multiomics-Arturo/methyl_data/metadata/ROSMAP_arrayMethylation_metaData.tsv")




# 1. AD data formatting (entrezIDs) for functional enrichment

# Methyl AD data
methyl_AD_genes <- array_meth_metadata %>%
                    filter(TargetID %in% rownames(stable_features_AD_data$methyl)) %>%
                    filter(!is.na(RefGene)) %>%
                    pull(RefGene) %>%
                    str_split(string = ., pattern = ";") %>% 
                    unlist

methyl_AD_entrezID <- bitr(
                        geneID   = methyl_AD_genes,
                        fromType = "SYMBOL", 
                        toType   = "ENTREZID", 
                        OrgDb    = org.Hs.eg.db
                      )
# 'select()' returned 1:1 mapping between keys and columns
# Aviso:
# In bitr(geneID = methyl_AD_genes, fromType = "SYMBOL", toType = "ENTREZID",  :
#  14.73% of input gene IDs are fail to map...

# Extract genes that failed to map
methyl_AD_genes_no_mapped <- setdiff(methyl_AD_genes,methyl_AD_entrezID$SYMBOL)
# length(methyl_AD_genes_no_mapped)
# [1] 450

# Get entrezID of genes that do not mapped with "ALIAS"
methyl_AD_entrezID_genes_no_mapped <- bitr(
                                        geneID   = methyl_AD_genes_no_mapped,
                                        fromType = "ALIAS", 
                                        toType   = "ENTREZID", 
                                        OrgDb    = org.Hs.eg.db
                                    )
# 'select()' returned 1:many mapping between keys and columns
# Aviso:
# In bitr(geneID = methyl_AD_genes_no_mapped, fromType = "ALIAS",  :
#  14.44% of input gene IDs are fail to map...
#
# length(setdiff(methyl_AD_genes_no_mapped, methyl_AD_entrezID_genes_no_mapped$ALIAS))
# [1] 65

# Remove duplicates
methyl_AD_entrezID_genes_no_mapped <- methyl_AD_entrezID_genes_no_mapped[!duplicated(methyl_AD_entrezID_genes_no_mapped$ALIAS),]

# Merge results
colnames(methyl_AD_entrezID_genes_no_mapped)[1] <- "SYMBOL"

methyl_AD_entrezID <- rbind(methyl_AD_entrezID, methyl_AD_entrezID_genes_no_mapped)
# dim(methyl_AD_entrezID)
# [1] 2991    2


# RNAseq AD data
rnaseq_AD_entrezID <- bitr(
                        geneID   = rownames(stable_features_AD_data$rnaseq), 
                        fromType = "ENSEMBL", 
                        toType   = "ENTREZID", 
                        OrgDb    = org.Hs.eg.db
                      )
# 'select()' returned 1:many mapping between keys and columns
# Aviso:
# In bitr(geneID = rownames(stable_features_AD_data$rnaseq), fromType = "ENSEMBL",  :
#  0.96% of input gene IDs are fail to map...

# setdiff(rownames(stable_features_AD_data$rnaseq), rnaseq_AD_entrezID$ENSEMBL)
# [1] "ENSG00000270607" "ENSG00000260816"


# miRNA AD data

mirna_AD_genes <- rownames(stable_features_AD_data$mirna)
# mirna_AD_genes
# [1] "hsa-miR-129-5p" "hsa-miR-132"    "hsa-miR-193b"   "hsa-miR-200a"
# [5] "hsa-miR-633"    "hsa-miR-744"    "hsa-miR-1260"   "hsa-miR-1275"

# HGNC format:
mirna_AD_genes <- c("MIR129-1", "MIR129-2", "MIR132", "MIR193B", "MIR200A", "MIR633", "MIR744", "MIR1260A", "MIR1275")

mirna_AD_entrezID <- bitr(
                        geneID   = mirna_AD_genes, 
                        fromType = "SYMBOL", 
                        toType   = "ENTREZID", 
                        OrgDb    = org.Hs.eg.db
                      )
# 'select()' returned 1:1 mapping between keys and columns

# Merge all entrezIDs (from each AD block)
AD_genes <- unique(c(methyl_AD_entrezID$ENTREZID, rnaseq_AD_entrezID$ENTREZID, mirna_AD_entrezID$ENTREZID))

length(AD_genes)
# [1] 3171 // For compareCluster



# 2. Control data formatting for functional enrichment

# Methyl control data
methyl_control_genes <- array_meth_metadata %>% 
                          filter(TargetID %in% rownames(stable_features_control_data$methyl)) %>%
                          filter(!is.na(RefGene)) %>%
                          pull(RefGene) %>%
                          str_split(string = ., pattern = ";") %>% 
                          unlist

methyl_control_entrezID <- bitr(
                            geneID   = methyl_control_genes,
                            fromType = "SYMBOL", 
                            toType   = "ENTREZID", 
                            OrgDb    = org.Hs.eg.db
                          )
# 'select()' returned 1:many mapping between keys and columns
# Aviso:
# In bitr(geneID = methyl_control_genes, fromType = "SYMBOL", toType = "ENTREZID",  :
#  14.59% of input gene IDs are fail to map...

# Remove duplicates
methyl_control_entrezID <- methyl_control_entrezID[!duplicated(methyl_control_entrezID$SYMBOL), ]

# Extract genes that failed to map
methyl_control_genes_no_mapped <- setdiff(methyl_control_genes, methyl_control_entrezID$SYMBOL)
# length(methyl_control_genes_no_mapped)
# [1] 468

# Get entrezID of genes that do not mapped with "ALIAS"
methyl_control_entrezID_genes_no_mapped <- bitr(
                                            geneID   = methyl_control_genes_no_mapped,
                                            fromType = "ALIAS", 
                                            toType   = "ENTREZID", 
                                            OrgDb    = org.Hs.eg.db
                                          )
# 'select()' returned 1:many mapping between keys and columns
# Aviso:
# In bitr(geneID = methyl_control_genes_no_mapped, fromType = "ALIAS",  :
#  14.74% of input gene IDs are fail to map...
#
# length(setdiff(methyl_control_genes_no_mapped, methyl_control_entrezID_genes_no_mapped$ALIAS))
# [1] 69

# Remove duplicates
methyl_control_entrezID_genes_no_mapped <- methyl_control_entrezID_genes_no_mapped[!duplicated(methyl_control_entrezID_genes_no_mapped$ALIAS), ]

# Merge results
colnames(methyl_control_entrezID_genes_no_mapped)[1] <- "SYMBOL"

methyl_control_entrezID <- rbind(methyl_control_entrezID, methyl_control_entrezID_genes_no_mapped)
# dim(methyl_control_entrezID)
# [1] 3139    2

# RNAseq control data
rnaseq_control_entrezID <- bitr(
                            geneID   = rownames(stable_features_control_data$rnaseq), 
                            fromType = "ENSEMBL", 
                            toType   = "ENTREZID", 
                            OrgDb    = org.Hs.eg.db
                          )
# 'select()' returned 1:many mapping between keys and columns
# Aviso:
# In bitr(geneID = rownames(stable_features_control_data$rnaseq),  :
#  4.37% of input gene IDs are fail to map...
#
# length(setdiff(rownames(stable_features_control_data$rnaseq), rnaseq_control_entrezID$ENSEMBL))
# [1] 9


# miRNA AD data

mirna_control_genes <- rownames(stable_features_control_data$mirna)
# mirna_control_genes
# [1] "hsa-miR-16"     "hsa-miR-103"    "hsa-miR-106b"   "hsa-miR-140-5p"
# [5] "hsa-miR-195"    "hsa-miR-497"    "hsa-miR-515-3p" "hsa-miR-627"
# [9] "ebv-miR-BART5"  "kshv-miR-K12-2"
# last 2 comes from viral sequences

# HGNC format:
mirna_control_genes <- c("MIR16-1", "MIR16-2", "MIR103A1","MIR103A2", "MIR106B", "MIR140", "MIR195", "MIR497", "MIR515-1", "MIR515-2", "MIR627")

mirna_control_entrezID <- bitr(
                            geneID   = mirna_control_genes, 
                            fromType = "SYMBOL", 
                            toType   = "ENTREZID", 
                            OrgDb    = org.Hs.eg.db
                          )
# 'select()' returned 1:1 mapping between keys and columns

# Merge all entrezIDs (from each control block)
control_genes <- unique(c(methyl_control_entrezID$ENTREZID, rnaseq_control_entrezID$ENTREZID, mirna_control_entrezID$ENTREZID))

length(control_genes)
# [1] 3318



# 3. Functional enrichment Analysis
enrich_BP <- compareCluster(
              geneClusters  = list(AD = AD_genes, Control = control_genes),
	            fun           = "enrichGO",
	            OrgDb         = org.Hs.eg.db,
	            ont           = "BP",
            	readable      = TRUE,
            	pAdjustMethod = "fdr",
              pvalueCutoff  = 0.05
          )

enrich_KEGG <- compareCluster(
              geneClusters  = list(AD = AD_genes, Control = control_genes),
	            fun           = "enrichKEGG",
	            organism      = "hsa",
            	pAdjustMethod = "fdr",
              pvalueCutoff  = 0.05
          )



enrich_BP_pval.01 <- compareCluster(
              geneClusters  = list(AD = AD_genes, Control = control_genes),
	            fun           = "enrichGO",
	            OrgDb         = org.Hs.eg.db,
	            ont           = "BP",
            	readable      = TRUE,
            	pAdjustMethod = "fdr",
              pvalueCutoff  = 0.01
          )

enrich_KEGG_pval.01 <- compareCluster(
              geneClusters  = list(AD = AD_genes, Control = control_genes),
	            fun           = "enrichKEGG",
	            organism      = "hsa",
            	pAdjustMethod = "fdr",
              pvalueCutoff  = 0.01
          )







