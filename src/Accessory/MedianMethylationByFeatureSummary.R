
## Script takes the count (coverage) information for all cpgs with at least 5x coverage starting with the methylKit object, subsets these
 # counts by feature and summarizes (takes the median) methylation for each individual by feature. The script stores the counts subset by 
 # feature as a nested list (count type with feature nested, i.e. tC$exon would give the total counts table for CpGs located within exons).
 # This script was also used to generate the median methylation by feature summary table with associated meta data for each individuals 
 # (i.e. treatment and time).

library(methylKit)
library(dplyr)
library(matrixStats)

## Path
inputDirCounts = "/shared_lab/20180226_RNAseq_2017OAExp/DNAm/processed_samples/03_CytoSummaries/methylKit_outputs"
inputDirFeature = "/shared_lab/20180226_RNAseq_2017OAExp/DNAm/processed_samples/09_CpG_summary"
outputFolder = "featureMethylationLevel"
saveDate="20200202"

# meta data
meta <- readRDS("/shared_lab/20180226_RNAseq_2017OAExp/DNAm/metadata/metadata_20190811.RData")
meta <- meta[meta$ID != "17099",]
# Read in Countsmetadata_20190811.RData
setwd(inputDirCounts)
refCoord <- readRDS("methylKitObj_all_cov5Filtered_united.RData")
coord <- paste0(refCoord$chr,"_",refCoord$start-1,"_",refCoord$end+1)
methCount <- read.csv("methylKit_cov5_counts/methylKitObj_all_cov5Filtered_united_MethylCCounts.csv")
totalCount <- read.csv("methylKit_cov5_counts/methylKitObj_all_cov5Filtered_united_totalCounts.csv")
methCount <- as.matrix(methCount)
totalCount <- as.matrix(totalCount)
class(methCount) <- "numeric"
class(totalCount) <- "numeric"
beta <- methCount/totalCount

methCount <- data.frame(as.numeric(methCount))
totalCount <- data.frame(as.numeric(totalCount))
beta <- data.frame(beta)

mC <- data.frame(coord,methCount)
colnames(mC) <- c("coord",meta$ID)
tC <- data.frame(coord,totalCount)
colnames(tC) <- c("coord",meta$ID)
b <- data.frame(coord,beta)
colnames(b) <- c("coord",meta$ID)

# Feature bed files
setwd(inputDirFeature)
exonF <- read.delim("20200129_CpG_cov5_Exon.txt",sep="\t",header=FALSE)
geneF <- read.delim("20200129_CpG_cov5_gene.txt",sep="\t",header=FALSE)
IntronF <- read.delim("20200129_CpG_cov5_Intron.txt",sep="\t",header=FALSE)
IntergenicF <- read.delim("20200129_CpG_cov5_Intergenic.txt",sep="\t",header=FALSE)

# Remove redundant CpGs in each feature (this can happen when exons overlap for instance)
exonF_unique <- exonF[!duplicated(paste0(exonF$V4,"_",exonF$V5,"_",exonF$V6)),4:6]
exonF_unique <- data.frame(coord=paste0(exonF_unique$V4,"_",exonF_unique$V5,"_",exonF_unique$V6))
geneF_unique <- geneF[!duplicated(paste0(geneF$V4,"_",geneF$V5,"_",geneF$V6)),4:6]
geneF_unique <- data.frame(coord=paste0(geneF_unique$V4,"_",geneF_unique$V5,"_",geneF_unique$V6))
IntronF_unique <- IntronF[!duplicated(paste0(IntronF$V4,"_",IntronF$V5,"_",IntronF$V6)),4:6]
IntronF_unique <- data.frame(coord=paste0(IntronF_unique$V4,"_",IntronF_unique$V5,"_",IntronF_unique$V6))
IntergenicF_unique <- IntergenicF[!duplicated(paste0(IntergenicF$V4,"_",IntergenicF$V5,"_",IntergenicF$V6)),4:6]
IntergenicF_unique <- data.frame(coord=paste0(IntergenicF_unique$V4,"_",IntergenicF_unique$V5,"_",IntergenicF_unique$V6))

# Overlap DNAm count and beta values with features

# Methylation Counts
exon_mC <- left_join(exonF_unique,mC) 
gene_mC <- left_join(geneF_unique,mC)
Intron_mC <- left_join(IntronF_unique,mC)
Intergenic_mC <- left_join(IntergenicF_unique,mC)
feature_mC <- list(exon=exon_mC,gene=gene_mC,Intron=Intron_mC,Intergenic=Intergenic_mC)
# Total Counts (coverage)
exon_tC <- left_join(exonF_unique,tC) 
gene_tC <- left_join(geneF_unique,tC)
Intron_tC <- left_join(IntronF_unique,tC)
Intergenic_tC <- left_join(IntergenicF_unique,tC)
feature_tC <- list(exon=exon_tC,gene=gene_tC,Intron=Intron_tC,Intergenic=Intergenic_tC)
# Beta
exon_beta <- left_join(exonF_unique,b) 
gene_beta <- left_join(geneF_unique,b)
Intron_beta <- left_join(IntronF_unique,b)
Intergenic_beta <- left_join(IntergenicF_unique,b)
feature_beta <- list(exon=exon_beta,gene=gene_beta,Intron=Intron_beta,Intergenic=Intergenic_beta)

feature_all <- list(mC=feature_mC,tC=feature_tC,beta=feature_beta)
# Summarize beta by treatment x time
# median of mean methylation among samples
med_meth <- NULL
ID_vec <- NULL
trt <- as.character(NULL)
time <- as.character(NULL)
feature <- NULL
featureNames <- c("Exon","gene","Intron","Intergenic")
for(i in 1:length(feature_beta)){
  med_meth <- c(med_meth,colMedians(as.matrix(feature_beta[[1]][,2:ncol(feature_beta[[1]])])))
  ID_vec <- c(ID_vec,meta$ID)
  trt <- c(trt,as.character(meta$Treatment))
  time <- c(time,as.character(meta$Time))
  feature <- c(feature,rep(featureNames[i],times=c(nrow(meta))))
}
medianSummary <- data.frame(feature=feature,treatment=trt,time=time,ID=ID_vec,median_methylation=med_meth)

# Save
if (file.exists(outputFolder)){
  setwd(file.path(inputDirFeature,outputFolder))
} else {
  dir.create(file.path(inputDirFeature,outputFolder))
  setwd(file.path(inputDirFeature,outputFolder))
}

saveRDS(feature_all,paste0(saveDate,"_AllCountsList_cov5_byFeature.RData"))
write.csv(medianSummary,paste0(saveDate,"medianMethylationByFeatureSummaryTable.csv"))

