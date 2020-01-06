#### Differential Methylation with bayesian mixed binomial model implemented using BRMS ####

## Packages
knitr::opts_chunk$set(echo = TRUE)
library(kableExtra)
library(dplyr)
library(reshape2)
library(stringr)
library(VennDiagram)
library(GSEABase)
library(RColorBrewer)
pal <- brewer.pal(n = 8, name = 'Dark2')

#### Data ####
# Sample Information
meta_c <- readRDS(paste0(wd,"meta/metadata_20190811.RData"))
meta <- meta_c[meta_c$ID != "17099",]
gene_ID <- readRDS(paste0(wd,"references/Exon_GeneLoc.RData"))
cds <- readRDS(paste0(wd,"references/CDS_wGoTerms_GeneLOC.RData"))
# Gene Expression 
gc <- readRDS("Transcriptomic/DGEListObj_withIndWeights_filterApproach2_plannedContrastMatrix.RData")
gc2 <- readRDS("Transcriptomic/gene_EBayesObj.RData")
gc3 <- readRDS("Transcriptomic/gene_postVoomAndNormalization_DGEListObj.RData")
gc4 <- readRDS("Transcriptomic/gene_preNormalization_DGEListObj.RData")
# Methylation Data
meth_models <- readRDS("DNAm/Final_DNAm_gene_BRMS_modelSummary.RData")
meth_marg <- readRDS("DNAm/Final_DNAm_gene_BRMS_modelMarginalEffects.RData")
d_meth <- readRDS("DNAm/Final_DNAm_gene_BRMS_plannedComparisons.RData")
d_meta <- readRDS("DNAm/CG_unstranded_summaryTable_geneOnly_5.RData")
d_beta <- readRDS("DNAm/CG_unstranded_beta_geneOnly_5.RData")
# Removing individual with low quality (located in day 8 trt ambient)
d_beta <-  d_beta[,colnames(d_beta) != "17099"]
LOC <- sub(".*Name=(.*?);.*","\\1",d_meta$attribute,perl=TRUE) # Extract LOC ID for each gene we have 
## Gene level summary information
cpgMean_all <-  read.table("DNAm/meanBetaPerFeature.txt",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")

#### Selecting out Significant Differentially Methylated Loci ####
## Using posterior probability to select for significant cpgs
alpha <- 20000
sum((1-d_meth$Post.Prob[,1]) < (0.05/alpha))
which((1-d_meth$Post.Prob[,1]) < (0.05/alpha))
sum((1-d_meth$Post.Prob[,2]) < (0.05/alpha))
which((1-d_meth$Post.Prob[,2]) < (0.05/alpha))
sum((1-d_meth$Post.Prob[,3]) < (0.05/alpha))
sum((1-d_meth$Post.Prob[,4]) < (0.05/alpha))
## Using Evid.Ratio
sum(log10(d_meth$Evid.Ratio[,1]) > 2)
er_index1 <- which(log10(d_meth$Evid.Ratio[,1]) > 2)
sum(log10(d_meth$Evid.Ratio[,2]) > 2)
er_index2 <- which(log10(d_meth$Evid.Ratio[,2]) > 2)
sum(log10(d_meth$Evid.Ratio[,3]) > 2)
er_index3 <- which(log10(d_meth$Evid.Ratio[,3]) > 2)
sum(log10(d_meth$Evid.Ratio[,4]) > 2)
er_index4 <- which(log10(d_meth$Evid.Ratio[,4]) > 2)
## Looking by credibility intervals
sum(d_meth$CI.Lower[,1] > 0)
ci_index1 <- which(d_meth$CI.Lower[,1] > 0)#trt
sum(d_meth$CI.Upper[,2] < 0)
ci_index2 <- which(d_meth$CI.Upper[,2] < 0)#trt
sum(d_meth$CI.Lower[,3] > 0)
ci_index3 <- which(d_meth$CI.Lower[,3] > 0)#time
sum(d_meth$CI.Upper[,4] < 0)
ci_index4 <- which(d_meth$CI.Upper[,4] < 0)#time
ci_index5 <- which(d_meth$CI.Lower[,5]>0 | d_meth$CI.Upper[,5]<0)#interaction

index1 <- er_index1[match(ci_index1,er_index1)]
index2 <- er_index2[match(ci_index2,er_index2)]
index_trt <- unique(index1,index2 )
index3 <- er_index3[match(ci_index3,er_index3)]
index4 <- er_index4[match(ci_index4,er_index4)]
index_time <- unique(index3,index4)

sig_design <- data.frame(rbind(cbind(index_trt,"Trt"),
                               cbind(index_time,"Time"),
                               cbind(ci_index5,"Interaction")))
colnames(sig_design) <- c("index","comp")
sig_design$index <- as.numeric(as.character(sig_design$index))
sig_design$comp <- as.character(sig_design$comp)
# This will transform long form data into matrix with n unique locations and rows equal to the number of primary comparisons (trt, time, interaction) = 1 is equal to significant for that variable.
primary_index <- dcast(sig_design,index~comp,length)
## Specific post hoc comparisons
#EvidRatio
colnames(d_meth$Evid.Ratio)
er_index1.1 <- which(log10(d_meth$Evid.Ratio[,6]) > 2)
er_index1.2 <- which(log10(d_meth$Evid.Ratio[,7]) > 2)
er_index1.3 <- which(log10(d_meth$Evid.Ratio[,8]) > 2)
er_index1.4 <- which(log10(d_meth$Evid.Ratio[,9]) > 2)
er_index1.5 <- which(log10(d_meth$Evid.Ratio[,10]) > 2)
er_index1.6 <- which(log10(d_meth$Evid.Ratio[,11]) > 2)
er_index1.7 <- which(log10(d_meth$Evid.Ratio[,12]) > 2)
er_index1.8 <- which(log10(d_meth$Evid.Ratio[,13]) > 2)
er_index1.9 <- which(log10(d_meth$Evid.Ratio[,14]) > 2)
er_index1.10 <- which(log10(d_meth$Evid.Ratio[,15]) > 2)
er_index1.11 <- which(log10(d_meth$Evid.Ratio[,16]) > 2)
er_index1.12 <- which(log10(d_meth$Evid.Ratio[,17]) > 2)
#Credibility Interval
colnames(d_meth$Estimate)
ci_index1.1 <- which(d_meth$CI.Lower[,6] > 0)#trtC9_E9
ci_index1.2 <- which(d_meth$CI.Upper[,7] < 0)#trtC9_E9
ci_index1.3 <- which(d_meth$CI.Lower[,8] > 0)#trtC9_C80
ci_index1.4 <- which(d_meth$CI.Upper[,9] < 0)#trtC9_C80
ci_index1.5 <- which(d_meth$CI.Lower[,10] > 0)#trtC9_E80
ci_index1.6 <- which(d_meth$CI.Upper[,11] < 0)#trtC9_E80
ci_index1.7 <- which(d_meth$CI.Lower[,12] > 0)#trtE9_C80
ci_index1.8 <- which(d_meth$CI.Upper[,13] < 0)#trtE9_C80
ci_index1.9 <- which(d_meth$CI.Lower[,14] > 0)#trtE9_E80
ci_index1.10 <- which(d_meth$CI.Upper[,15] < 0)#trtE9_E80
ci_index1.11<- which(d_meth$CI.Lower[,16] > 0)#trtC80_E80
ci_index1.12 <- which(d_meth$CI.Upper[,17] < 0)#trtC80_E80

index1.1 <- Reduce(intersect, list(ci_index1.1,er_index1.1))
index1.2 <- Reduce(intersect,list(ci_index1.2,er_index1.2)) 
index1.3 <- Reduce(intersect,list(ci_index1.3,er_index1.3))
index1.4 <- Reduce(intersect,list(ci_index1.4,er_index1.4))
index1.5 <- Reduce(intersect,list(ci_index1.5,er_index1.5))
index1.6 <- Reduce(intersect,list(ci_index1.6,er_index1.6))
index1.7 <- Reduce(intersect,list(ci_index1.7,er_index1.7))
index1.8 <- Reduce(intersect,list(ci_index1.8,er_index1.8))
index1.9 <- Reduce(intersect,list(ci_index1.9,er_index1.9))
index1.10 <- Reduce(intersect,list(ci_index1.10,er_index1.10))
index1.11 <- Reduce(intersect,list(ci_index1.11,er_index1.11))
index1.12 <- Reduce(intersect,list(ci_index1.12,er_index1.12))

trtC9_E9 <- unique(index1.1,index1.2)
trtC9_C80 <- unique(index1.3,index1.4)
trtC9_E80 <- unique(index1.5,index1.6)
trtE9_C80 <- unique(index1.7,index1.8)
trtE9_E80 <- unique(index1.9,index1.10)
trtC80_E80 <- unique(index1.11,index1.12)

sig_design2 <- data.frame(rbind(
  cbind(ci_index5,"Interaction"),
  cbind(trtC9_E9,"trtC9_E9"),
  cbind(trtC9_C80,"trtC9_C80"),
  cbind(trtC9_E80,"trtC9_E80"),
  cbind(trtE9_C80,"trtE9_C80"),
  cbind(trtE9_E80,"trtE9_E80"),
  cbind(trtC80_E80,"trtC80_E80")))

colnames(sig_design2) <- c("index","comp")
sig_design2$index <- as.numeric(as.character(sig_design2$index))
sig_design2$comp <- as.character(sig_design2$comp)

#### Supp Figure - Venn Diagram ####
ID <-  d_meta$ID
loci_names <- data.frame(ID=ID,index=1:length(ID))
primary_index_save <- left_join(primary_index,loci_names)

sig_pos_trt <- primary_index$index[primary_index$Trt==1]
sig_pos_time <- primary_index$index[primary_index$Time==1]
sig_pos_inter <- primary_index$index[primary_index$Interaction==1]

Complete_overlap <- Reduce(intersect, list(sig_pos_trt,sig_pos_time,sig_pos_inter))
trt_time_overlap <- Reduce(intersect, list(sig_pos_trt,sig_pos_time))
trt_inter_overlap <- Reduce(intersect, list(sig_pos_trt,sig_pos_inter))
time_inter_overlap <- Reduce(intersect, list(sig_pos_time,sig_pos_inter))

grid.newpage()
draw.triple.venn(area1 = length(sig_pos_trt),
                 area2 = length(sig_pos_time),
                 area3 = length(sig_pos_inter),
                 n12 = 1, #length(trt_time_overlap),
                 n23 = length(time_inter_overlap),
                 n13 = 1, #length(trt_inter_overlap),
                 n123 = 1, #length(Complete_overlap),
                 category = c("Treatment", "Time", "Interaction"),
                 fill = c("white","white","white"),
                 lty = c(1,2,3),
                 alpha = 0.6,
                 cex = 2.5,
                 cat.dist = c(0.07,0.06,0.04),
                 cat.cex = c(2.8,2.8,2.8))

#### Selecting out Loci Significant by trt (main effect) ####

## Extract information about the significant treatment genes
trt_indexes <- primary_index$index[primary_index$Trt==1] 
gene_labels <- gc$genes[match(LOC[trt_indexes],gc$genes$GENEID),]
# Calculate log fold change (using raw counts)
gene_Expression <- gc4$counts[match(LOC[trt_indexes],gc$genes$GENEID),]
gene_Expression_sum_C <- rowMeans(gene_Expression[,meta_c$treatment == 400])
gene_Expression_sum_E <- rowMeans(gene_Expression[,meta_c$treatment == 2800])
logFoldChange <- log(gene_Expression_sum_E/gene_Expression_sum_C)
# Coding sequence labels
cds_labels <- cds[match(LOC[trt_indexes],cds$gene_id),]
# CpG Meta Data
trt_meta<-d_meta[trt_indexes,] 
# CpG % change
trt_C_meanBeta <- rowMeans(d_beta[trt_indexes,meta$treatment==2800])
trt_E_meanBeta <- rowMeans(d_beta[trt_indexes,meta$treatment==400])
finalPercentChange <- round((trt_E_meanBeta-trt_C_meanBeta)*100,2)

#write.csv(cbind(gene_labels,cds_labels,logFoldChange,
#                trt_meta,trt_C_meanBeta,trt_E_meanBeta,
#                finalPercentChange),
#          "DNAm/Final_diffMethlyation_Treatment_raw.csv")

#### Significant out Loci Significant by trt (post hoc at either time point) #### 
ph_sig <- data.frame(rbind(
  cbind(ci_index5,"Interaction"),
  cbind(index1.1,"C9_E9_Hypo"),
  cbind(index1.2,"C9_E9_Hyper"),
  cbind(index1.11,"C80_E80_Hypo"),
  cbind(index1.12,"C80_E80_Hyper")))

colnames(ph_sig) <- c("index","comp")
ph_sig$index <- as.numeric(as.character(ph_sig$index))
ph_sig$comp <- as.character(ph_sig$comp)

ph_index_trt <- dcast(ph_sig,index~comp,length)
ph_index_trt <- ph_index_trt[ph_index_trt$Interaction==1,]

ph_f <- ph_index_trt[rowSums(ph_index_trt[,2:5])>0,]
ph_f <- ph_f[,-6]

ph_f_09_index <- rowSums(ph_f[,4:5])>0
ph_f_80_index <- rowSums(ph_f[,2:3])>0

## Summarize methylation by treatment for significant loci
# Day 9
d_beta_trt <- d_beta[ph_f$index,]
d_beta_trt_09_C <- rowMeans(d_beta_trt[ph_f_09_index,meta$SFV == "09.400"])
d_beta_trt_09_E <- rowMeans(d_beta_trt[ph_f_09_index,meta$SFV == "09.2800"])
d_beta_trt_09_diff <- c(d_beta_trt_09_E-d_beta_trt_09_C)*100
d_beta_trt_09_sum <- data.frame(Timepoint="Day_9",
                                index = ph_f$index[ph_f_09_index],
                                Mean_Ambient = as.vector(d_beta_trt_09_C),
                                Mean_OA = as.vector(d_beta_trt_09_E),
                                Change = as.vector(d_beta_trt_09_diff))
# Day 80
d_beta_trt_80_C <- rowMeans(d_beta_trt[ph_f_80_index,meta$SFV == "80.400"])
d_beta_trt_80_E <- rowMeans(d_beta_trt[ph_f_80_index,meta$SFV == "80.2800"])
d_beta_trt_80_diff <- c(d_beta_trt_80_E-d_beta_trt_80_C)*100
d_beta_trt_80_sum <- data.frame(Timepoint="Day_80",
                                index = ph_f$index[ph_f_80_index],
                                Mean_Ambient = as.vector(d_beta_trt_80_C),
                                Mean_OA = as.vector(d_beta_trt_80_E),
                                Change = as.vector(d_beta_trt_80_diff))
d_beta_trt_sum <- rbind(d_beta_trt_09_sum,d_beta_trt_80_sum)
#gene_labels_f <- gc$genes[match(LOC[d_beta_trt_sum$index],gc$genes$GENEID),]
cds_labels_f <- cds[match(LOC[d_beta_trt_sum$index],cds$gene_id),]
col <- c("protein_id","gene_id","Description","GO.SLIM","InterproScan.GO.ID")
cds_labels_f_red <- subset(cds_labels_f,select=col)
comb <- cbind(d_beta_trt_sum,cds_labels_f_red)
trt_meta_ph<-d_meta[comb$index,]
comb <- cbind(comb,trt_meta_ph)
## Gene Expression for corresponding genes
# Day 9
#gene_Expression_9 <- gc$E[match(LOC[ph_f$index[ph_f_09_index]],gc$genes$GENEID),]
gene_Expression_9 <- gc4$counts[match(LOC[ph_f$index[ph_f_09_index]],gc$genes$GENEID),]
gene_Expression_9[is.na(gene_Expression_9[,1]),] <- 0
gene_Expression_sum_C9 <- rowMeans(gene_Expression_9[,meta_c$SFV == "09.400"])
gene_Expression_sum_E9 <- rowMeans(gene_Expression_9[,meta_c$SFV == "09.2800"])
logFoldChange_9 <- log(gene_Expression_sum_E9/gene_Expression_sum_C9)
# Day 80 
#gene_Expression_80<- gc$E[match(LOC[ph_f$index[ph_f_80_index]],gc$genes$GENEID),]
gene_Expression_80<- gc4$counts[match(LOC[ph_f$index[ph_f_80_index]],gc$genes$GENEID),]
gene_Expression_80[is.na(gene_Expression_80[,1]),] <- 0
gene_Expression_sum_C80 <- rowMeans(gene_Expression_80[,meta_c$SFV == "80.400"])
gene_Expression_sum_E80 <- rowMeans(gene_Expression_80[,meta_c$SFV == "80.2800"])
logFoldChange_80 <- log(gene_Expression_sum_E80/gene_Expression_sum_C80)
# Combine the two timepoints together
logFoldChange <- c(logFoldChange_9,logFoldChange_80)

#Combine DNA MEthylation and GE logFoldChange
comb$logFoldChange <- logFoldChange

## GO slim term categories ##
# Separate out the go terms
split_go <- strsplit(comb$GO.SLIM,split = ";")
go_list <- matrix(ncol=2,nrow=length(unlist(split_go)))
go_list <- data.frame(go_list)
colnames(go_list) <- c("index","go_slim_raw")
m <- 1
for(i in 1:length(split_go)){
  if(length(split_go[[i]])>0){
    temp <- split_go[[i]]
    for(j in 1:length(temp)){
      go_list[m,1] <- comb$index[i]
      go_list[m,2] <- temp[j]
      m <- c(m+1)
    }
  }
}

go_list$go_slim_final <- str_sub(go_list$go_slim_raw,start = -10,end =-1)
# Lets consider the most basic go term for each cpgs
go_list_first <- go_list[!duplicated(go_list$index), ]

# NOTE: this may take time to run
#Get biological function
myCollection <- GOCollection(go_list$go_slim_final)
#fl <- system.file("extdata", "goslim_plant.obo", package="GSEABase")
fl <- "transcriptomic/go.obo"
slim <- getOBOCollection(fl)
goSlim_output_MF <- goSlim(myCollection, slim, "MF")
goSlim_output_BP <- goSlim(myCollection, slim, "BP")
goSlim_output_CC <- goSlim(myCollection, slim, "CC")
go_MF <- data.frame(go_slim_final=rownames(goSlim_output_MF),MF_term=goSlim_output_MF$Term)
go_BP <- data.frame(go_slim_final=rownames(goSlim_output_BP),BP_term=goSlim_output_BP$Term)
go_CC <- data.frame(go_slim_final=rownames(goSlim_output_CC),CC_term=goSlim_output_CC$Term)

go_list_mf <- left_join(go_list,go_MF)
go_list_mf <- go_list_mf[!is.na(go_list_mf$MF_term),]
go_list_bp <- left_join(go_list,go_BP)
go_list_bp <- go_list_bp[!is.na(go_list_bp$BP_term),]
go_list_cc <- left_join(go_list,go_CC)
go_list_cc <- go_list_cc[!is.na(go_list_cc$CC_term),]

length(ph_f$index)
go_list_mf_order <- go_list_mf[match(ph_f$index,go_list_mf$index),c(2:4)]
ph_f_update <- cbind(ph_f,go_list_mf_order)
go_list_bp_order <- go_list_bp[match(ph_f$index,go_list_bp$index),c(2:4)]
ph_f_update <- cbind(ph_f_update,go_list_bp_order)
go_list_cc_order <- go_list_cc[match(ph_f$index,go_list_cc$index),c(2:4)]
ph_f_update <- cbind(ph_f_update,go_list_cc_order)
# Creating a single variable for time
ph_f_update$diffMethylationTime <- "Empty"
ph_f_update$diffMethylationTime[ph_f_update$C80_E80_Hyper == 1 | ph_f_update$C80_E80_Hypo == 1] <- "Day_80"
ph_f_update$diffMethylationTime[ph_f_update$C9_E9_Hyper == 1 | ph_f_update$C9_E9_Hypo == 1] <- "Day_9"
# Creating a single variable for direction of change (hyper vs hypo)
ph_f_update$directionChange <- "Empty"
ph_f_update$directionChange[ph_f_update$C9_E9_Hyper == 1 | ph_f_update$C80_E80_Hyper == 1] <- "Hyper"
ph_f_update$directionChange[ph_f_update$C80_E80_Hypo == 1 | ph_f_update$C9_E9_Hypo == 1] <- "Hypo"
#Combining time and direction into a single variable
ph_f_update$singleMethylVariable <- paste0(ph_f_update$diffMethylationTime,"_",ph_f_update$directionChange)

ph_f_update_MF <- ph_f_update[!is.na(ph_f_update$MF_term),]
ph_f_update_BP <- ph_f_update[!is.na(ph_f_update$BP_term),]
ph_f_update_CC <- ph_f_update[!is.na(ph_f_update$CC_term),]

mf_table <- as.matrix(table(as.character(ph_f_update_MF$MF_term),ph_f_update_MF$singleMethylVariable))
bp_table <- as.matrix(table(as.character(ph_f_update_BP$BP_term),ph_f_update_BP$singleMethylVariable))
cc_table <- as.matrix(table(as.character(ph_f_update_CC$CC_term),ph_f_update_CC$singleMethylVariable))

## Order samples to combine them into single dataframe
comb_order <- comb[order(comb$index),]  
ph_f_update_order <- ph_f_update[order(ph_f_update$index),]
final_ph <- cbind(comb_order,ph_f_update_order)
#write.csv(final_ph,"DNAm/Final_diffMethlyation_PostHoc_raw.csv")