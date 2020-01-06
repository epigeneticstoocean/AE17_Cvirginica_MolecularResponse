#### Differential Gene Expresison ####

## packages
library(matrixStats,quietly = TRUE)
library(edgeR,quietly = TRUE)
library(limma,quietly = TRUE)
library(fdrtool,quietly = TRUE)

#### Data ####
wd <- "~/AE17_Cvirginica_MolecularResponse/data/Analysis"
# This should be set to the path for the local version of the `2017OAExp_Oysters` github repo.
readRDS(paste0(wd,"gene_postVoomAndNormalization_DGEListObj.RData"))

#### Analysis ####

#### Identify correlation between factors in design contrasts with blocking factor ####
gene_a2_o1_cor <- duplicateCorrelation(dge_gene_a2_o1_voom, design, block = meta$tankID)

#### Fitting Model ####
lmf_gene_a2_o1 <- lmFit(dge_gene_a2_o1_voom, design,
                        block = meta$tankID,
                        correlation = gene_a2_o1_cor$consensus.correlation)
#Refitting Option 1 with contrasts
lmf_gene_a2_o1_wContr <- contrasts.fit(lmf_gene_a2_o1,contr_mat)
##Run empiricial bayes protocol
bayes_gene_a2_o1_contr <- eBayes(lmf_gene_a2_o1_wContr,robust=TRUE)

## Saving Files
#saveRDS(bayes_gene_a2_o1_contr,paste0(wd,"/gene_EBayesObj.RData"))
