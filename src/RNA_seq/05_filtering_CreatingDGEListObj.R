#### Filtering and Storing Expression Matrices ####

## Packages 
library(edgeR,quietly = TRUE)
library(limma,quietly = TRUE)
#library(devtools)
#install_github("vqv/ggbiplot")
library(ggbiplot,quietly = TRUE)

#### Data#####
## Set working directory
wd <- "~/Github/AE17_Cvirginica_MolecularResponse/data"
outputDir <- "~/Github/AE17_Cvirginica_MolecularResponse/data/Analysis"
# This should be set to the path for the local version of the `2017OAExp_Oysters` github repo.
## RSEM counts 
# Gene matrix
RSEM <-  readRDS(paste0(wd,"/RNAseq/RSEM_gene_Summary.Rdata"))
# Separate out RSEM counts and rename rows with LOC ID
rsem_c <- RSEM$Count # Stick with raw gene counts  (not TPM or some other normalizaed count)
rm(RSEM)
## Transcript File  
# Transcript file
tran <- readRDS(paste0(wd,"/references/STAR_gnomon_tximportGeneFile.RData"))
# Genes in File
gene <- tran[!duplicated(tran$GENEID),]
## Meta Data  
# metaTemp <- read.csv("~/Github/AE17_Cvirginica_MolecularResponse/data/meta/AE17_RNAmetaData.csv",stringsAsFactors = FALSE)
# metaTemp$Treatment <- as.factor(metaTemp$Treatment)
# metaTemp$Time <- as.factor(metaTemp$Time)
# metaTemp$SFV <- as.factor(metaTemp$SFV)
# metaTemp$Pop <- as.factor(metaTemp$Pop)
#saveRDS(metaTemp,"~/Github/AE17_Cvirginica_MolecularResponse/data/meta/AE17_RNAmetaData.RData")
meta <- readRDS(paste0(wd,"/meta/AE17_RNAmetaData.RData"))
#Create new factor levels (one for each level combination)
meta$SFVrn <- as.factor(paste0("D",meta$SFV))
meta$Sample_Index <- as.factor(meta$sample_index)
meta$TankID <- as.factor(meta$tankID)

#### Data Manipulation ####
# Order genes from annotation list to match the order in the count matrix
gene_order <- gene[order(gene$gene_id),]
identical(rownames(rsem_c),gene_order$gene_id) # TRUE confirms the order matches
# Relabel the rows of the count matrix with LOC ID
rownames(rsem_c) <- gene_order$GENEID

# Reorder by LOCID 
geneC_full <- rsem_c

rm(rsem_c)

#### Filtering counts ####
## Round to whole counts
geneC_all <- round(geneC_full)

## Genes 
# Breaking down expression coverage by treatment*time combination
#Day 9 Trt 2800
keep_D9.2800 <- rowSums(cpm(geneC_all[,meta$SFVrn=="D09.2800"])>=1) >= 5
sum(keep_D9.2800)
#Day 9 Trt 400
keep_D9.400 <- rowSums(cpm(geneC_all[,meta$SFVrn=="D09.400"])>=1) >= 5
sum(keep_D9.400)
#Day 80 Trt 2800
keep_D80.2800 <- rowSums(cpm(geneC_all[,meta$SFVrn=="D80.2800"])>=1) >= 5
sum(keep_D80.2800)
#Day 80 Trt 400
keep_D80.400 <- rowSums(cpm(geneC_all[,meta$SFVrn=="D80.400"])>=1) >= 5
sum(keep_D80.400)

keep_gene_a2 <- rowSums(cbind(keep_D9.2800,keep_D9.400,
                              keep_D80.2800,keep_D80.400)) >= 1

# Filter 
geneC_a2 <- geneC_all[keep_gene_a2, ]
gene_final_a2 <- gene_order[keep_gene_a2,]

## Create DGEList
dge_gene_a2 <- DGEList(geneC_a2,genes = gene_final_a2) # counts - rsem

## Count summary tables
# All genes captures before filter
perSampleMeanCount <- mean(colSums(geneC_all)) #mean
perSampleSDCount <- sd(colSums(geneC_all)) #sd
perSampleMinCount <- min(colSums(geneC_all)) #min
perSampleMaxCount <- max(colSums(geneC_all)) #max
numberOfGenes <- sum(rowSums(geneC_all)>0) # Number of genes in dataset (rows)
table <- data.frame(Data="All counts",mean=perSampleMeanCount,sd=perSampleSDCount,min=perSampleMinCount,max=perSampleMaxCount,genes=numberOfGenes)
# Genes counts after filtering
perSampleMeanCountF <- mean(colSums(geneC_a2)) #mean
perSampleSDCountF <- sd(colSums(geneC_a2)) #sd
perSampleMinCountF <- min(colSums(geneC_a2)) #min
perSampleMaxCountF <- max(colSums(geneC_a2)) #max
numberOfGenesF <- sum(rowSums(geneC_a2)>0) # Number of genes in dataset (rows)
table <- rbind(table,
               data.frame(Data="Filtered counts",
                          mean=perSampleMeanCountF,sd=perSampleSDCountF,
                          min=perSampleMinCountF,max=perSampleMaxCountF,
                          genes=numberOfGenesF))
print(table)

#### Save initial DGEList objects####
#saveRDS(dge_gene_a2,paste0(outputDir,"/RNA_gene_preNormalization_DGEListObj.RData"))

#### Normalization with edgeR ####
# Calculate normalization factors for scaling raw lib. size
dge_gene_a2_norm <- calcNormFactors(dge_gene_a2,method = "TMMwsp") # gene - approach 2
# Bar plot of normalization factors
barplot(dge_gene_a2_norm$samples$norm.factors~rownames(dge_gene_a2_norm$samples),
        las=2,ylab="Normalization factor",xlab="Samples")
plotMDS(dge_gene_a2_norm, col = as.numeric(meta$SFVrn))

## Create design matrix  
design <- model.matrix(~0+SFVrn,data=meta) # 0+ is needed here otherwise the first level defaults to 1.
#Rename columns
colnames(design) <- levels(meta$SFVrn)

#### Transform and create observational level weights ####
## Gene Features 
dge_gene_a2_o1_voom <- voomWithQualityWeights(dge_gene_a2_norm,design,plot = TRUE)
## Plots
barplot(dge_gene_a2_o1_voom$targets$sample.weights~rownames(dge_gene_a2_o1_voom$targets),
        las=2,ylab="Sample Specific Weights",xlab="Samples")
ge.pca <- prcomp(t(dge_gene_a2_o1_voom$E), center = TRUE,scale = TRUE)
ggbiplot(ge.pca,
         ellipse=FALSE,
         obs.scale = 1,
         var.scale = 1,
         var.axes=FALSE,
         labels= colnames(dge_gene_a2_o1_voom$E),
         groups=meta$SFVrn)
plotMDS(dge_gene_a2_o1_voom, col = as.numeric(meta$SFVrn))

# Saving the final transformation of the data (after individual weights)
# saveRDS(dge_gene_a2_o1_voom,paste0(outputDir,"/RNA_gene_postVoomAndNormalization_DGEListObj.RData"))

# Dendrogram based on co-expression clustering. Further illustrates 17005 as an outlier.
out <- dist(t(scale(dge_gene_a2_o1_voom$E)))
out2 <- hclust(out,method = "complete")
plot(out2)