#### Script creates RData objects of the phenotype data for downstream analysis ####

#### Data ####
setwd("/home/downeyam/Github/AE17_Cvirginica_MolecularResponse/data/Phenotype/")
# EPF data
pheno <- readRDS("AE17_summaryPhenotype_alltimepoints.RData") 
# Final version of bouyant weights
cal <- read.csv("AE17_BouyantWeight_Final.csv",stringsAsFactors = FALSE)
cal$ID <- as.integer(cal$ID) # missing coercion message due to oysters that were already dead
# Final version of sample data
samp <- read.csv("AE17_CollectionInfo.csv")

# Remove oysters that died
samp <- samp[samp$Acc_Exp != "DEAD",]
# Only consider oysters from the exposure (remove data from acclimation and additional point)
samp <- samp[samp$Acc_Exp == "Exposure",]
samp <- samp[samp$tpnum_actual < 10,]
ID <- data.frame(ID=samp$ID)

# Subset EPF and calcification data to include only oysters from exposure
pheno <- left_join(ID,pheno)
cal <- left_join(ID,cal)
cal <- left_join(ID,cal)

### Save EPF timeseries data (used for Fig1
saveRDS(pheno,"AE17_summaryPhenotype_exposure.RData")
out <- read.csv("AE17_summaryPhenotype_exposure.csv")
out <- out[,-1]
colNames <- data.frame(num=c(1:ncol(out)),cols=colnames(out))
#### Long term EPF pH and calcification data ####
colnames(pheno)
columns <- c("ID","pCO2_fac","Tank_pH","EPF_pH","Daily_Pcnt_Change_Corrected","EPF_envAdj")
pheno_target <- subset(pheno,select=columns)
# Merge with sample sheet to include site info
samp_red <- subset(samp,select=c(1,9,38,39,40,41,42))
cal <- subset(cal,select=c(1:14,15,17,21,23,27,29,39:51))
cal <- merge(samp_red,cal)
  
#### Calcification Data ####
setwd("/home/downeyam/Github/AE17_Cvirginica_MolecularResponse/data/Phenotype/")
# Final version of bouyant weights
cal<-read.csv("AE17_BouyantWeight_Final.csv")
# Final version of sample data
samp<-read.csv("AE17_CollectionInfo.csv")
# Final EPF fluid data 
pheno <- readRDS("AE17_summaryPhenotype_exposure.RData")
pheno_red <- pheno[pheno$timepoint != 81,]
pheno_red <- pheno[pheno$timepoint > 40,]
pheno_red <- subset(pheno_red,select=c(1,11,21,30,37,44))

# Merge with sample sheet to include site info
samp_red <- subset(samp,select=c(1,9,38,39,40,41,42))
cal <- subset(cal,select=c(1:14,15,17,21,23,27,29,39:51))
cal <- merge(samp_red,cal)

# Filter complete dataset
# Remove individuals that don't have values for the third bouyant weight
cal_red <- cal[!is.na(cal$bw2_bw3_daily_sizeCorr),]
# Remove the extra oysters used for the calcein experiment and from tp 81
cal_red2 <-  cal_red[cal_red$sample_date < 20170824,]

## Calculate calcification
cal_red2$calcification <- cal_red2$bw2_bw3_daily/cal_red2$dry2_final