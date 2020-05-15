# Short script is used to combine the CpG-By-Gene counts and summary statistics into a single data.frame
library(dplyr)
library(matrixStats)

## Pathes
inputDir <- "~/Github/AE17_Cvirginica_MolecularResponse/"
outputDir <- "~/Github/AE17_Cvirginica_MolecularResponse/data/MBDBS_seq/20200130_CpGbyGeneSummary/"
# Should be set to folder with all CpGByGene Summary txt files
setwd(inputDir)

# Gene level summary count for all CpGs in the genome
aC <- read.delim("data/MBDBS_seq/20200130_CpGbyGeneSummary/20200130_allCpG_Count.txt",sep="\t",header=FALSE)
aC_rev <- data.frame(label=paste0(aC$V1,"_",aC$V2,"_",aC$V3),chr=aC$V1,start=aC$V2,end=aC$V3,all_count=aC$V4)

# Gene level summary count for all CpGs with at least 5x for each individuals
covC <- read.delim("data/MBDBS_seq/20200130_CpGbyGeneSummary/20200130_cov5CpG_Count.txt",sep="\t",header=FALSE)
covC_rev <- data.frame(label=paste0(covC$V1,"_",covC$V2,"_",covC$V3),cov5_count=covC$V4)
# Gene level summary stats for all CpGs with at least 5x for each individuals
summaryC <- read.delim("data/MBDBS_seq/20200130_CpGbyGeneSummary/20200130_cov5CpG_Summarize.txt",sep="\t",header=FALSE)
summaryC_rev <- data.frame(label=paste0(summaryC$V1,"_",summaryC$V2,"_",summaryC$V3),summaryC[,4:13])
summaryC_label <-c("label","Mean","Mean_9C","Mean_9E","Mean_80C","Mean_80E","cv_9C","cv_mean_AmongTrt","diff_Trt","diff_tp9_Trt","diff_tp80_Trt")
colnames(summaryC_rev) <- summaryC_label

# Gene level summary count for all diff. methylated CpGs by OA treatment (all timepoints)
DML <- read.delim("data/MBDBS_seq/20200130_CpGbyGeneSummary/20200130_DMLTrt_Count.txt",sep="\t",header=FALSE)
DML_rev <- data.frame(label=paste0(DML$V1,"_",DML$V2,"_",DML$V3),DML_count=DML$V4)

# Gene level summary count for diff. methylated CpGs by OA treatment (Day 9)
DML_tp9 <- read.delim("data/MBDBS_seq/20200130_CpGbyGeneSummary/20200130_DML_tp9_Count.txt",sep="\t",header=FALSE)
DML_tp9_rev <- data.frame(label=paste0(DML_tp9$V1,"_",DML_tp9$V2,"_",DML_tp9$V3),DML_tp9_count=DML_tp9$V4)

# Gene level summary count for diff. methylated CpGs by OA treatment (Day 80)
DML_tp80 <- read.delim("data/MBDBS_seq/20200130_CpGbyGeneSummary/20200130_DML_tp80_Count.txt",sep="\t",header=FALSE)
DML_tp80_rev <- data.frame(label=paste0(DML_tp80$V1,"_",DML_tp80$V2,"_",DML_tp80$V3),DML_tp80_count=DML_tp80$V4)

# Series of left joines to add the CpG count and summary information
comb <- left_join(aC_rev,covC_rev)
comb <- left_join(comb,DML_rev)
comb <- left_join(comb,DML_tp9_rev)
comb <- left_join(comb,DML_tp80_rev)
comb <- left_join(comb,summaryC_rev)
# Save combined data.frame in an .RData format (can be read into R with readRDS()) in the same folder at the original CpG
# count and summary txt files.
saveRDS(comb,paste0(outputDir,"gene_CpGcoverageSummary.RData"))

