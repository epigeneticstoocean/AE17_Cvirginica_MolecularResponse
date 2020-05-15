
## Libraries
library(edgeR)
library(matrixStats)

#### Data ####

## path
inputDir <- "~/Github/AE17_Cvirginica_MolecularResponse/"
# set directory
setwd(inputDir)

# Read in the gene level CpG count and summary data.frame
cpgS <- readRDS("data/MBDBS_seq/20200130_CpGbyGeneSummary/gene_CpGcoverageSummary.RData")

# Read in the gene expression data
ts<- readRDS("data/Analysis/RNA_gene_preNormalization_DGEListObj.RData")

e <- ts$counts
e_cpm_log <- log2(cpm(e)+ 1) # log2 cpm transformation
e <- e_cpm_log[,colnames(e_cpm_log) != "17099"] # remove the 17099 here because CpG is missing this ind.

# Number of exons per gene
exons <- read.delim("data/MBDBS_seq/20200130_CpGbyGeneSummary/20200130_exonInGene_Count.txt",sep="\t",header=FALSE)
exon_labels <- paste0(exons$V1,"_",exons$V2,"_",exons$V3)
# Exon counts and gene lengths
exons_rev <- data.frame(label=exon_labels,exon_count=exons$V4,gene_length=exons$V3-exons$V2)

# Read in the gene reference file - this contains the chromosome and starting gene coordinate (used to 
# uniquely identify genes in the CpG file) and the gene_id (LOCXXXXX) which will be used to merge the expression
# and CpG data.
ref <- readRDS("data/references/gene_GeneLoc.RData")
meta <- readRDS("data/meta/metadata_20190811.RData")
meta_min <- meta[meta$ID != "17099",]

#### Summarize the gene expression data ####
## Create a summary of methylation for each CpG ##
# Mean gene expression all sample3s
mean <- rowMeans(e)
# Mean gene expression Control day 9
mean_9C <- rowMeans(e[,meta_min$SFV == "09.400"])
# Mean gene expression Exposed day 9
mean_9E <- rowMeans(e[,meta_min$SFV == "09.2800"])
# Mean gene expression Control day 80
mean_80C <- rowMeans(e[,meta_min$SFV == "80.400"])
# Mean gene expression Exposed day 80
mean_80E <- rowMeans(e[,meta_min$SFV == "80.2800"])
# Coefficient of Variation for Control day 9 samples
sd_9C <- rowSds(e[,meta_min$SFV == "09.400"])
cv_9C <- sd_9C/mean_9C
# Coefficient of Variation for Exposed day 9 samples
sd_9E <- rowSds(e[,meta_min$SFV == "09.2800"])
cv_9E <- sd_9E/mean_9E
# Coefficient of Variation for Control day 80 samples
sd_80C <- rowSds(e[,meta_min$SFV == "80.400"])
cv_80C <- sd_9C/mean_9C
# Coefficient of Variation for Exposed day 80 samples
sd_80E <- rowSds(e[,meta_min$SFV == "80.2800"])
cv_80E <- sd_80E/mean_80E
# Coefficient of Variation among treatments
meanAmongTrt <- cbind(mean_9C,mean_9E,mean_80C,mean_80E)
sd_meanAmongTrt <- rowSds(meanAmongTrt)
mean_meanAmongTrt <- rowMeans(meanAmongTrt)
cv_mean_AmongTrt <- sd_meanAmongTrt/mean_meanAmongTrt
# Difference in OA Trt
mean_C <- rowMeans(e[,meta_min$Treatment ==400])
mean_E <- rowMeans(e[,meta_min$Treatment == 2800])
diff_Trt <- mean_E-mean_C
# Difference in gene expression tp9
diff_tp9_Trt <- mean_9E-mean_9C
# Difference in gene expression tp80
diff_tp80_Trt <- mean_80E-mean_80C
# Creat matrix of all CpG summary stats
gene_summary <- cbind(mean,mean_9C,mean_9E,mean_80C,mean_80E,
                      cv_9C,cv_9E,cv_80C,cv_80E,cv_mean_AmongTrt,
                      diff_Trt,diff_tp9_Trt,diff_tp80_Trt)
gene_summary[is.na(gene_summary)]<-0
colnames(gene_summary) <- paste0("gene_",colnames(gene_summary))
gene_summary <- data.frame(gene_id = rownames(gene_summary),gene_summary)
gene_annot <- left_join(gene_summary,ref)
gene_annot$label <- paste0(gene_annot$chr,"_",gene_annot$start,"_",gene_annot$end)
gene_reduce <- data.frame(label=gene_annot$label,gene_annot[,2:14])


# Combine gene expression, exon count and CpG (DNA mehtylation data)
gSum <- left_join(cpgS,gene_reduce)
gSum <- left_join(gSum,exons_rev)
# Reducing the data.frame to include places were we have all the data.
gSum_red<-gSum[!is.na(gSum$cov5_count),]
gSum_red<-gSum_red[!is.na(gSum_red$gene_mean),]
# Saving reduced data.frame
saveRDS(gSum,"data/Analysis/Multi_geneSummaryComplete.RData")
# Saving the reduced data.frame
saveRDS(gSum_red,"data/Analysis/Multi_geneSummaryReduced.RData")
