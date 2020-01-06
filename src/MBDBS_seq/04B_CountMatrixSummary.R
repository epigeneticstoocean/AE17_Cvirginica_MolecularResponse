#### Libraries ####
library(dplyr)
library(ggplot2)
library(matrixStats)
library(gridExtra)
library(R.utils)

sprintf("Reading in data")
#### Read in Data ####
gene <- readRDS("/shared_lab/20180226_RNAseq_2017OAExp/DNAm/reference/gene_GeneLoc.RData")
exon <- readRDS("/shared_lab/20180226_RNAseq_2017OAExp/DNAm/reference/Exon_GeneLoc.RData")

setwd("/shared_lab/20180226_RNAseq_2017OAExp/DNAm/trimmedSamples/singleSample_trimScript_20190719/countMatrix")
# All cytosines 
meta <-readRDS("All_CytoSum_summaryTable.RData")
meta$index <- c(1:nrow(meta))
# All cytosines within genes
meta_gene <- readRDS("CytoSummary_Gene_Match.RData")
## TEMPORARY TAKE OUT WHEN RUNNING REAL (make sure none have index value larger than the number of rows in subset data)
#meta_gene <- meta_gene[mex < 10000,]
meta_gene_m1 <- meta_gene[meta_gene$match_num==1,]
# Some matching coding issues in an earlier step create a few none matching rows rather than completely
# omitting them, here we just remove incomplete cases (should only be a few rows.
meta_gene_m1_noNA <- meta_gene_m1[which(!is.na(meta_gene_m1$cg_index)),]
meta_gene_m1_noNA <- meta_gene_m1[which(!is.na(meta_gene_m1$transcript_id)),]
meta_gene <- meta_gene_m1_noNA

# All cytosines within exons
meta_exon <- readRDS("CytoSummary_Exon_Match.RData")
#meta_exon <- meta_exon[meta_exon$cg_index < 10000,]
# Need to only look at the first match for exons since each exon can show up
# more than once if there are mutiple isoforms. This can be seen when mutiple Parent
# rnaIDs exist for a single gene, and those different rnas (isoforms) share an exon.
# This will show up redundantly in the the meta_exon df, so just selecting out the first
# match will select only the first usage of the exon, which for exon level analysis is fine
meta_exon_m1 <- meta_exon[meta_exon$match_num==1,]
# Some matching coding issues in an earlier step create a few none matching rows rather than completely 
# omitting them, here we just remove incomplete cases (should only be a few rows.
meta_exon_m1_noNA <- meta_exon_m1[which(!is.na(meta_exon_m1$cg_index)),]
meta_exon_m1_noNA <- meta_exon_m1[which(!is.na(meta_exon_m1$transcript_id)),]
meta_exon <- meta_exon_m1_noNA

## Count matrices 
# Methylated
mC <- readRDS("All_CytoSum_methylCountMatrix.RData")
# Unmethylated
umC <- readRDS("All_CytoSum_unmethylCountMatrix.RData")
# Total
tC <- readRDS("All_CytoSum_TotalCountMatrix.RData")
## Beta (proportion methylated)
beta <-readRDS("All_CytoSum_BetaMatrix.RData")
beta[is.na(beta)] <- 0 

## Removing sample 17099 because of poor sequencing (this should increase the retained number of loci 
tC  <- tC[colnames(tC)!=17099,]
mC  <- mC[colnames(tC)!=17099,]
umC <- umC[colnames(tC)!=17099,]
beta <- beta[colnames(tC)!=17099,]

# Directory for outputs
#setwd("/shared_lab/20180226_RNAseq_2017OAExp/DNAm/trimmedSamples/singleSample_trimScript_20190719/countSummary")
setwd("/shared_lab/20180226_RNAseq_2017OAExp/DNAm/trimmedSamples/singleSample_trimScript_20190719/countSummary_V2")



sprintf("Pruning by min count threshold...")
# All cytosine
# Pruning by variable min count thresholds
tC_min <- apply(tC,1,min) # vector of minimum counts for each cpg
tC_median <- apply(tC,1,median) # vector of median counts for each CpG

# Summing each CpG
tC_sum <- rowSums(tC)
## Sum 5 threshold that Yaamini used
meta_sum5 <- meta[tC_sum>=5,]
mC_sum5 <- mC[tC_sum>=5,]
umC_sum5 <- umC[tC_sum>=5,]
tC_sum5 <- tC[tC_sum>=5,]
beta_sum5 <- beta[tC_sum>=5,]
saveRDS(meta_sum5,"Final_meta_sum5.RData")
saveRDS(mC_sum5,"Final_mC_sum5.RData")
saveRDS(umC_sum5,"Final_umC_sum5.RData")
saveRDS(tC_sum5,"Final_tC_sum5.RData")
saveRDS(beta_sum5,"Final_beta_sum5.RData")

### Minimum threshold sample
#At least 1 count per sample
meta_0 <- meta[tC_min>0,]
mC_0 <- mC[tC_min>0,]
umC_0 <- umC[tC_min>0,]
tC_0 <- tC[tC_min>0,]
beta_0 <- beta[tC_min>0,]

# 5 counts per sample
meta_5 <- meta[tC_min>=5,]
mC_5 <- mC[tC_min>=5,]
umC_5 <- umC[tC_min>=5,]
tC_5 <- tC[tC_min>=5,]
beta_5 <- beta[tC_min>=5,]

# 10 counts per sample
meta_10 <- meta[tC_min>=10,]
mC_10 <- mC[tC_min>=10,]
umC_10 <- umC[tC_min>=10,]
tC_10 <- tC[tC_min>=10,]
beta_10 <- beta[tC_min>=10,]

# Save outputs
saveRDS(meta_5,"Final_meta_5.RData")
saveRDS(mC_5,"Final_mC_5.RData")
saveRDS(umC_5,"Final_umC_5.RData")
saveRDS(tC_5,"Final_tC_5.RData")
saveRDS(beta_5,"Final_beta_5.RData")
saveRDS(meta_10,"Final_meta_10.RData")
saveRDS(mC_10,"Final_mC_10.RData")
saveRDS(umC_10,"Final_umC_10.RData")
saveRDS(tC_10,"Final_tC_10.RData")
saveRDS(beta_10,"Final_beta_10.RData")

sprintf("Matching with genes....")
# Cytosines within genes - variable thresholds

mC_gene <- mC[meta_gene$cg_index,]
umC_gene <- umC[meta_gene$cg_index,]
tC_gene <- tC[meta_gene$cg_index,]
beta_gene <- beta[meta_gene$cg_index,]
# Prune by variabel min count thresholds
tC_gene_min <- apply(tC_gene,1,min) # vector of minimum counts for each cpg
# At least 1 per sample
meta_gene_0 <- meta_gene[tC_gene_min>0,]
mC_gene_0 <- mC_gene[tC_gene_min>0,]
umC_gene_0 <- umC_gene[tC_gene_min>0,]
tC_gene_0 <- tC_gene[tC_gene_min>0,]
beta_gene_0 <- beta_gene[tC_gene_min>0,]
# At least 5 per sample
meta_gene_5 <- meta_gene[tC_gene_min>=5,]
mC_gene_5 <- mC_gene[tC_gene_min>=5,]
umC_gene_5 <- umC_gene[tC_gene_min>=5,]
tC_gene_5 <- tC_gene[tC_gene_min>=5,]
beta_gene_5 <- beta_gene[tC_gene_min>=5,]
# At least 10 per sample
meta_gene_10 <- meta_gene[tC_gene_min>=10,]
mC_gene_10 <- mC_gene[tC_gene_min>=10,]
umC_gene_10 <- umC_gene[tC_gene_min>=10,]
tC_gene_10 <- tC_gene[tC_gene_min>=10,]
beta_gene_10 <- beta_gene[tC_gene_min>=10,]

saveRDS(meta_gene_5,"Final_meta_gene_5.RData")
saveRDS(mC_gene_5,"Final_mC_gene_5.RData")
saveRDS(umC_gene_5,"Final_umC_gene_5.RData")
saveRDS(tC_gene_5,"Final_tC_gene_5.RData")
saveRDS(beta_gene_5,"Final_beta_gene_5.RData")
saveRDS(meta_gene_10,"Final_meta_gene_10.RData")
saveRDS(mC_gene_10,"Final_mC_gene_10.RData")
saveRDS(umC_gene_10,"Final_umC_gene_10.RData")
saveRDS(tC_gene_10,"Final_tC_gene_10.RData")
saveRDS(beta_gene_10,"Final_beta_gene_10.RData")

sprintf("Matching with exons ...")
# Cytosines with Exon - variable thresholds

mC_exon <- mC[meta_exon$cg_index,]
umC_exon <- umC[meta_exon$cg_index,]
tC_exon <- tC[meta_exon$cg_index,]
beta_exon <- beta[meta_exon$cg_index,]
# Prune by variabel min count thresholds
tC_exon_min <- apply(tC_exon,1,min) # vector of minimum counts for each cpg
# At least 1 per sample
meta_exon_0 <- meta_exon[tC_exon_min>0,]
mC_exon_0 <- mC_exon[tC_exon_min>0,]
umC_exon_0 <- umC_exon[tC_exon_min>0,]
tC_exon_0 <- tC_exon[tC_exon_min>0,]
beta_exon_0 <- beta_exon[tC_exon_min>0,]
# At least 5 per sample
meta_exon_5 <- meta_exon[tC_exon_min>=5,]
mC_exon_5 <- mC_exon[tC_exon_min>=5,]
umC_exon_5 <- umC_exon[tC_exon_min>=5,]
tC_exon_5 <- tC_exon[tC_exon_min>=5,]
beta_exon_5 <- beta_exon[tC_exon_min>=5,]
# At least 10 per sample
meta_exon_10 <- meta_exon[tC_exon_min>=10,]
mC_exon_10 <- mC_exon[tC_exon_min>=10,]
umC_exon_10 <- umC_exon[tC_exon_min>=10,]
tC_exon_10 <- tC_exon[tC_exon_min>=10,]
beta_exon_10 <- beta_exon[tC_exon_min>=10,]

saveRDS(meta_exon_5,"Final_meta_exon_5.RData")
saveRDS(mC_exon_5,"Final_mC_exon_5.RData")
saveRDS(umC_exon_5,"Final_umC_exon_5.RData")
saveRDS(tC_exon_5,"Final_tC_exon_5.RData")
saveRDS(beta_exon_5,"Final_beta_exon_5.RData")
saveRDS(meta_exon_10,"Final_meta_exon_10.RData")
saveRDS(mC_exon_10,"Final_mC_exon_10.RData")
saveRDS(umC_exon_10,"Final_umC_exon_10.RData")
saveRDS(tC_exon_10,"Final_tC_exon_10.RData")
saveRDS(beta_exon_10,"Final_beta_exon_10.RData")

## NOTE: This reads in a count matrix from the bedgraph option
 # which by default only considers CpGs

sprintf("Creating summaries...")
#### Summaries ####

# Descriptive Methylation 
# All CpGs
beta_mean <- rowMeans(beta)
beta_sd <- unlist(array(beta,1,sd))
# All CpGs with counts
beta_0_mean <- rowMeans(beta_0)
beta_0_sd <- unlist(array(beta_0,1,sd))
beta_0_CV <- beta_0_sd/beta_0_mean
beta_5_mean <- rowMeans(beta_5)
beta_5_sd <- unlist(array(beta_5,1,sd))
beta_5_CV <- beta_5_sd/beta_5_mean
beta_10_mean <- rowMeans(beta_10)
beta_10_sd <- unlist(array(beta_10,1,sd))
beta_10_CV <- beta_10_sd/beta_10_mean
### CpGs in genes
beta_gene_0_mean <- rowMeans(beta_gene_0)
beta_gene_0_sd <- unlist(array(beta_gene_0,1,sd))
beta_gene_0_CV <- beta_gene_0_sd/beta_gene_0_mean
beta_gene_5_mean <- rowMeans(beta_gene_5)
beta_gene_5_sd <- unlist(array(beta_gene_5,1,sd))
beta_gene_5_CV <- beta_gene_5_sd/beta_gene_5_mean
beta_gene_10_mean <- rowMeans(beta_gene_10)
beta_gene_10_sd <- unlist(array(beta_gene_10,1,sd))
beta_gene_10_CV <- beta_gene_10_sd/beta_gene_10_mean
## CpG w/ 10threshold - split by condition and time
beta_gene_10_control_9 <- beta_gene_10[,c(3,5,6,17,18,19)]
beta_gene_10_control_80 <- beta_gene_10[,c(1,2,8,11,22,23)]
beta_gene_10_exposed_9 <- beta_gene_10[,c(7,9,12,13,15,21)]
beta_gene_10_exposed_80 <- beta_gene_10[,c(4,10,14,16,20,24)]
# CpG w/ 10 threshold sumary split by condition and time 
# Control - Day 9
beta_gene_10_control_9_mean <- rowMeans(beta_gene_10_control_9)
beta_gene_10_control_9_sd <- unlist(array(beta_gene_10_control_9,1,sd))
beta_gene_10_control_9_CV <- beta_gene_10_control_9_sd/beta_gene_10_control_9_mean
# Control - Day 80
beta_gene_10_control_80_mean <- rowMeans(beta_gene_10_control_80)
beta_gene_10_control_80_sd <- unlist(array(beta_gene_10_control_80,1,sd))
beta_gene_10_control_80_CV <- beta_gene_10_control_80_sd/beta_gene_10_control_80_mean
# Exposed - Day 9
beta_gene_10_exposed_9_mean <- rowMeans(beta_gene_10_exposed_9)
beta_gene_10_exposed_9_sd <- unlist(array(beta_gene_10_exposed_9,1,sd))
beta_gene_10_exposed_9_CV <- beta_gene_10_exposed_9_sd/beta_gene_10_exposed_9_mean
# Exposed - Day 80
beta_gene_10_exposed_80_mean <- rowMeans(beta_gene_10_exposed_80)
beta_gene_10_exposed_80_sd <- unlist(array(beta_gene_10_exposed_80,1,sd))
beta_gene_10_exposed_80_CV <- beta_gene_10_exposed_80_sd/beta_gene_10_exposed_80_mean

# Gene Summary Table
gene_summary <- cbind(beta_gene_10_mean,beta_gene_10_sd,beta_gene_10_CV,
                      beta_gene_10_control_9_mean,beta_gene_10_control_9_sd,beta_gene_10_control_9_CV,
                      beta_gene_10_control_80_mean,beta_gene_10_control_80_sd,beta_gene_10_control_80_CV,
                      beta_gene_10_exposed_9_mean,beta_gene_10_exposed_9_sd,beta_gene_10_exposed_9_CV,
                      beta_gene_10_exposed_80_mean,beta_gene_10_exposed_80_sd,beta_gene_10_exposed_80_CV)
saveRDS(gene_summary,"Final_beta_gene_summary.RData")

### CpGs in exons
beta_exon_0_mean <- rowMeans(beta_exon_0)
beta_exon_0_sd <- unlist(array(beta_exon_0,1,sd))
beta_exon_0_CV <- beta_exon_0_sd/beta_exon_0_mean
beta_exon_5_mean <- rowMeans(beta_exon_5)
beta_exon_5_sd <- unlist(array(beta_exon_5,1,sd))
beta_exon_5_CV <- beta_exon_5_sd/beta_exon_5_mean
beta_exon_10_mean <- rowMeans(beta_exon_10)
beta_exon_10_sd <- unlist(array(beta_exon_10,1,sd))
beta_exon_10_CV <- beta_exon_10_sd/beta_exon_10_mean
## CpG w/ 10threshold - split by condition and time
beta_exon_10_control_9 <- beta_exon_10[,c(3,5,6,17,18,19)]
beta_exon_10_control_80 <- beta_exon_10[,c(1,2,8,11,22,23)]
beta_exon_10_exposed_9 <- beta_exon_10[,c(7,9,12,13,15,21)]
beta_exon_10_exposed_80 <- beta_exon_10[,c(4,10,14,16,20,24)]
# CpG w/ 10 threshold sumary split by condition and time 
# Control - Day 9
beta_exon_10_control_9_mean <- rowMeans(beta_exon_10_control_9)
beta_exon_10_control_9_sd <- unlist(array(beta_exon_10_control_9,1,sd))
beta_exon_10_control_9_CV <- beta_exon_10_control_9_sd/beta_exon_10_control_9_mean
# Control - Day 80
beta_exon_10_control_80_mean <- rowMeans(beta_exon_10_control_80)
beta_exon_10_control_80_sd <- unlist(array(beta_exon_10_control_80,1,sd))
beta_exon_10_control_80_CV <- beta_exon_10_control_80_sd/beta_exon_10_control_80_mean
# Exposed - Day 9
beta_exon_10_exposed_9_mean <- rowMeans(beta_exon_10_exposed_9)
beta_exon_10_exposed_9_sd <- unlist(array(beta_exon_10_exposed_9,1,sd))
beta_exon_10_exposed_9_CV <- beta_exon_10_exposed_9_sd/beta_exon_10_exposed_9_mean
# Exposed - Day 80
beta_exon_10_exposed_80_mean <- rowMeans(beta_exon_10_exposed_80)
beta_exon_10_exposed_80_sd <- unlist(array(beta_exon_10_exposed_80,1,sd))
beta_exon_10_exposed_80_CV <- beta_exon_10_exposed_80_sd/beta_exon_10_exposed_80_mean

# Gene Summary Table
exon_summary <- cbind(beta_exon_10_mean,beta_exon_10_sd,beta_exon_10_CV,
                      beta_exon_10_control_9_mean,beta_exon_10_control_9_sd,beta_exon_10_control_9_CV,
                      beta_exon_10_control_80_mean,beta_exon_10_control_80_sd,beta_exon_10_control_80_CV,
                      beta_exon_10_exposed_9_mean,beta_exon_10_exposed_9_sd,beta_exon_10_exposed_9_CV,
                      beta_exon_10_exposed_80_mean,beta_exon_10_exposed_80_sd,beta_exon_10_exposed_80_CV)
saveRDS(exon_summary,"Final_beta_exon_summary.RData")

sprintf("Creating summary table...")
#### Summary Table ####

nam <- c("Total_CpG",
         "CpG_0","CpG_5x","CpG_10x",
         "CpG_5x_gene","CpG_10x_gene",
         "CpG_5x_exon","CpG_10x_exon")

counts<-c(nrow(tC),
          nrow(tC_0),nrow(tC_5),nrow(tC_10),
          nrow(tC_gene_5),nrow(tC_gene_10),
          nrow(tC_exon_5),nrow(tC_exon_10))

mean_prop_methylation <- c(0.0,
                           mean(beta_0_mean),mean(beta_5_mean),mean(beta_10_mean),
                           mean(beta_gene_5_mean),mean(beta_gene_10_mean),
                           mean(beta_exon_5_mean),mean(beta_exon_10_mean))
sd_prop_methylation <- c(0.0,
                         sd(beta_0_mean),sd(beta_5_mean),sd(beta_10_mean),
                         sd(beta_gene_5_mean),sd(beta_gene_10_mean),
                         sd(beta_exon_5_mean),sd(beta_exon_10_mean))

#png("Full_summaryTable.png")
counts_table <- cbind(nam,counts,mean_prop_methylation,sd_prop_methylation)
saveRDS(counts_table,"Simple_Methyl_Prop_Summary.RData")
#p <- kableExtra::kable(counts_table) %>% kableExtra::kable_styling()
#grid.arrange(p)
#dev.off()

sprintf("Printing figures...")
#### Figures ####
## Density plots for different thresholds
# With all CpGs
comb_plot_density <- data.frame(rbind(cbind(beta=beta_mean,Threshold="All"),
                           cbind(beta=beta_0_mean,Threshold="All : Coverage > 0"),
                           cbind(beta=beta_5_mean,Threshold="All : Coverage >= 5"),
                           cbind(beta=beta_10_mean,Threshold="All : Coverage >= 10"),
                           cbind(beta=beta_gene_0_mean,Threshold="Genes : Coverage > 0"),
                           cbind(beta=beta_gene_5_mean,Threshold="Genes : Coverage >= 5"),
                           cbind(beta=beta_gene_10_mean,Threshold="Genes : Coverage >= 10")
))
comb_plot_density$beta <- as.numeric(as.character(comb_plot_density$beta))

ggplot(comb_plot_density,aes(beta,colour=Threshold)) + 
  geom_density() + labs(title = "Mean Methylation Prop")
ggsave("Mean_Methyl_Prop_VariousThresholds.png")


# With just genes
#Means
comb_genes_Mean <- data.frame(rbind(cbind(beta=beta_gene_0_mean,Threshold="Genes : Coverage > 0"),
                                      cbind(beta=beta_gene_5_mean,Threshold="Genes : Coverage >= 5"),
                                      cbind(beta=beta_gene_10_mean,Threshold="Genes : Coverage >= 10")
))
comb_genes_Mean$beta <- as.numeric(as.character(comb_genes_Mean$beta))

ggplot(comb_genes_Mean,aes(beta,colour=Threshold)) + 
  geom_density() + xlim(0,1) +
  labs(title = "Mean Methylation Prop for just genes",xlab="Mean Methylation") 
ggsave("Mean_Methyl_Prop_Genes.png")

#SD
comb_genes_SD <- data.frame(rbind(cbind(beta=beta_gene_0_sd,Threshold="Genes : Coverage > 0"),
                                    cbind(beta=beta_gene_5_sd,Threshold="Genes : Coverage >= 5"),
                                    cbind(beta=beta_gene_10_sd,Threshold="Genes : Coverage >= 10")
))
comb_genes_SD$beta <- as.numeric(as.character(comb_genes_SD$beta))
ggplot(comb_genes_SD,aes(beta,colour=Threshold)) + 
  geom_density() + labs(title = "SD Methylation Prop for just genes") + xlim(0,1)
ggsave("SD_Mean_Methyl_Prop_Genes.png")

## Genes by treatment x time combination
comb_genesByTreatment_Mean <- data.frame(rbind(cbind(beta=beta_gene_10_control_9_mean,Threshold="Control : 9"),
                                  cbind(beta=beta_gene_10_control_80_mean,Threshold="Control : 80"),
                                  cbind(beta=beta_gene_10_exposed_9_mean,Threshold="Exposed : 9"),
                                  cbind(beta=beta_gene_10_exposed_80_mean,Threshold="Exposed : 80"))
)
comb_genesByTreatment_Mean$beta <- as.numeric(as.character(comb_genesByTreatment_Mean$beta))

ggplot(comb_genesByTreatment_Mean,aes(beta,colour=Threshold)) + 
  geom_density() + xlim(0,1) +
  labs(title = "Mean Methylation Prop for just genes by treatment x time (threshold >= 10)",xlab="Mean Methylation") 
ggsave("Mean_Methyl_Prop_Genes_DiffTreatments.png")

## Bar plot of genic methylation levels

comb_GlobalMean <- data.frame(rbind(cbind(beta=mean(beta_gene_10_control_9_mean),sd=sd(beta_gene_10_control_9_mean),Threshold="Control : 09"),
                                    cbind(beta=mean(beta_gene_10_control_80_mean),sd=sd(beta_gene_10_control_80_mean),Threshold="Control : 80"),
                                    cbind(beta=mean(beta_gene_10_exposed_9_mean),sd=sd(beta_gene_10_exposed_9_mean),Threshold="Exposed : 09"),
                                    cbind(beta=mean(beta_gene_10_exposed_80_mean),sd=sd(beta_gene_10_exposed_80_mean),Threshold="Exposed : 80"),
                                    cbind(beta=mean(beta_gene_10_mean),sd=sd(beta_gene_10_mean),Threshold=" All Genes "),
                                    cbind(beta=mean(beta_10_mean),sd=sd(beta_10_mean),Threshold="All CpGs ")
                                               ))
  comb_GlobalMean$beta <- as.numeric(as.character(comb_GlobalMean$beta))
comb_GlobalMean$sd <- as.numeric(as.character(comb_GlobalMean$sd))
ggplot(comb_GlobalMean,aes(x=Threshold,y=beta)) + 
  geom_errorbar(aes(ymin=beta-sd, ymax=beta+sd), width=.1) +
  geom_point() + 
  labs(title = "Global Mean Methylation Prop (threshold >= 10) ") + 
  ylim(0,1)
ggsave("Mean_Methyl_Prop_GlobalMean.png")


