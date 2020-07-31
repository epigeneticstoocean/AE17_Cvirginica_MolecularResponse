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
# Isoform (transcript) matrix
RSEM_t <-  readRDS(paste0(wd,"/RNAseq/RSEM_isoform_Summary.Rdata"))
rsem_t_c <- RSEM_t$Count
rm(RSEM_t)
## Transcript File  
# Transcript file
tran <- readRDS(paste0(wd,"/references/STAR_gnomon_tximportGeneFile.RData"))
# Genes in File
gene <- tran[!duplicated(tran$GENEID),]
## Meta Data  
meta <- readRDS(paste0(wd,"/meta/metadata_20190811.RData"))
meta$sampl_nameSimple <- substr(meta$sample_name,start = 4,stop=9)
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

# Order transcripts from annotation list to match the order from the count matrix
rsem_t_order <- rsem_t_c[order(rownames(rsem_t_c)),]
tran_order <- tran[order(tran$TXNAME),]
identical(rownames(rsem_t_order),tran_order$TXNAME)
# I am not renaming the rows here since the TXNAME will be unique, while the gene LOC will be redundant for isoforms from the same gene. However, the reordering here is still needed when the annotations get added to the DGEList later.

# Reorder by LOCID 
geneC_full <- rsem_c
tranC_full <- rsem_t_order

rm(rsem_c)
rm(rsem_t_order)

#### Filtering counts ####
## Round to whole counts
geneC_all <- round(geneC_full)
tranC_all <- round(tranC_full)

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

#### Not in the manuscript ####

## Script below was used to first remove outliers before processing which genes to filter and 
 # calculating sample corrections. This turns out not to impact the final interpretation (i.e.
 # that there are no differentially expressed genes), so it was removed from the supplement.

#### Outlier Removal ##
# List of samples removed
sample_list <- c("17005") #Consistent top outlier removed from data after normalization
# Alter meta file
meta_alt <- meta[!(meta$ID %in% sample_list),]
## Round to whole counts
geneC_alt <- geneC_all[,!(colnames(geneC_all) %in% sample_list)]
tranC_alt <- tranC_all[,!(colnames(tranC_all) %in% sample_list)]

## Genes 
# Breaking down expression coverage by treatment*time combination
#Day 9 Trt 2800

keep_D9.2800 <- rowSums(cpm(geneC_alt[,meta_alt$SFVrn=="D09.2800"])>=1) >= c(sum(meta_alt$SFVrn=="D09.2800")-1)
sum(keep_D9.2800)
#Day 9 Trt 400
keep_D9.400 <- rowSums(cpm(geneC_alt[,meta_alt$SFVrn=="D09.400"])>=1) >= c(sum(meta_alt$SFVrn=="D09.400")-1)
sum(keep_D9.400)
#Day 80 Trt 2800
keep_D80.2800 <- rowSums(cpm(geneC_alt[,meta_alt$SFVrn=="D80.2800"])>=1) >= c(sum(meta_alt$SFVrn=="D80.2800")-1)
sum(keep_D80.2800)
#Day 80 Trt 400
keep_D80.400 <- rowSums(cpm(geneC_alt[,meta_alt$SFVrn=="D80.400"])>=1) >= c(sum(meta_alt$SFVrn=="D80.400")-1)
sum(keep_D80.400)

keep_gene_a2 <- rowSums(cbind(keep_D9.2800,keep_D9.400,
                              keep_D80.2800,keep_D80.400)) >= 1

# Filter 
geneC_a2 <- geneC_alt[keep_gene_a2, ]
gene_final_a2 <- gene_order[keep_gene_a2,]

## Create DGEList
dge_gene_a2 <- DGEList(geneC_a2,genes = gene_final_a2) # counts - rsem

## Count summary tables
# All genes captures before filter
perSampleMeanCount <- mean(colSums(geneC_alt)) #mean
perSampleSDCount <- sd(colSums(geneC_alt)) #sd
perSampleMinCount <- min(colSums(geneC_alt)) #min
perSampleMaxCount <- max(colSums(geneC_alt)) #max
numberOfGenes <- sum(rowSums(geneC_alt)>0) # Number of genes in dataset (rows)
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

#### Save initial DGEList objects##
#saveRDS(dge_gene_a2,paste0(outputDir,"/RNA_gene_preNormalization_DGEListObj.RData"))

#### Normalization with edgeR ##
# Calculate normalization factors for scaling raw lib. size
dge_gene_a2_norm <- calcNormFactors(dge_gene_a2,method = "TMMwsp") # gene - approach 2
# Bar plot of normalization factors
barplot(dge_gene_a2_norm$samples$norm.factors~rownames(dge_gene_a2_norm$samples),
        las=2,ylab="Normalization factor",xlab="Samples")

## Create design matrix  
design <- model.matrix(~0+SFVrn,data=meta_alt) # 0+ is needed here otherwise the first level defaults to 1.
#Rename columns
colnames(design) <- levels(meta_alt$SFVrn)

#### Transform and create observational level weights ##
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
         groups=meta_alt$SFVrn)
# Saving the final transformation of the data (after individual weights)
saveRDS(dge_gene_a2_o1_voom,paste0(outputDir,"/RNA_gene_postVoomAndNormalization_DGEListObj_remove17005.RData"))
