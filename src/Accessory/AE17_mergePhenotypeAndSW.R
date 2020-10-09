#### Script for combining raw phenotype (EPF) and WC ####

# Created by: Alan Downey-Wall

# Purpose: Takes raw phenotype and water chemistry files and comebines them for downstream analysis.
#   Specifically, it:
#(i)  Summarizes water chemistry in three ways. First, is summarized tank chemistry for each individ  ual
#     over the entire duration of the exposure (net WC). Second, it summarizes the previous 2 weeks WC
#     for each individual. Third, for the first 33 days of the exposure (for samples collected on day 50 or 80).
#     The second summary (2 weeks) was used to calculate the delta pH (EPF pH - seawater pH). The third 
#     summary was used for looking at the relationship between pCO2 and calcification rate (estimated at day 33).
#(ii) Combines all samples points into `CompletePhenotype_AllSamples_final2020.csv` file. This
#     file contains individuals from both the AE17 exposure and a separate calcein and EPF recharge exp.
#     conducted by Louise cameron. Moreover, it also contains individuals that died during the exp.
#(iii)Cleaned data to to create `` file.
#     This file contains only individuals from the AE17 experiment, specifically those that survived and were
#     included in the exposure (i.e. we removed dead individuals, individuals from acclimation sampling,
#     and individuals from the extra sampling at the end of the experiment). A timepoint column was added based 
#     on sampling data.

#### Library and working directory ####
library(dplyr)
setwd("/home/downeyam/Github/AE17_Cvirginica_MolecularResponse")
#### Calcification Data ###
# Calcification Rates calculated from Buoyant Wgt data
cal <- read.delim("data/Phenotype/AE17_CalcificationInfo_ExperimentalExposureSamples_Simple.csv",sep=",")
#### Extra-pallial fluid pH data ###
# EPF pH data for all individuals using micro pH probe
epf <- read.delim("data/Phenotype/AE17_EPFpHData_ExperimentalExposureSamples_Simple.csv",sep=",")
#### Complete Chemistry Carb. Chemistry of the EPF ###
# Complete carbonate chemistry for a sub-sample of individuals from experiment
# Samples should largely correspond with individuals used for molecular work (day 9 and 80)
# Complete carb chem was calculated by Louise Cameron using measured DIC and pH.
#

#### Sample Data ###
# Oyster sample data for all individuals (even dead oysters and those from calcein exp.)
samples <- read.delim("data/Phenotype/AE17_metadata_ExperimentalExposureSamples_Simple.csv",sep=",")
samples <- samples[,-2]
#### Water Chemistry - Weekly ###
# Triweekly measurements of temperature, salinity, and pH (currently NBS scale)
wc_w <- read.delim("data/water_chem/AE17_WeeklySW_ExperimentalExposureSamples_Simple.csv",
                   sep=",",stringsAsFactors = FALSE)
#### Water Chemistry - Full Carb Chemistry.
# Complete carb. chem water chem. based on bottle samples run on VINDTA and carb. chem calculated
# using DIC and TA (when possible, pH substituted on occasion) and CO2Sys.
wc_c <- read.delim("data/water_chem/AE17_CarbChemSW_ExperimentalExposureSamples_Simple.csv",sep=",",
                   stringsAsFactors = FALSE)

#### Pruning DATA ####
# Removing columns from calcification .csv and epf carb. chemistry that overlap with other files.

## Calcification spreadsheet
# Columns removed: 
# "Rank","WaterSample_ID","shelf","tank","pCO2","DIC","TA","pCO2_calc","Calcite",
# "pH_meas","survival","survival_date","sample_date"
cal_sel <- cal[,c(2,9:10)]

# EPF complete carb chem spreadsheet
# Columns removed:
# "shelf" "tank"  "pCO2"
epf_sel <- epf[,c(1,7:10)]

#### Summarizing Water chemistry ####

# Create sampling date vector (needed for summarizing WC)
sample_date <- paste0(substring(samples$sample_date,1,4),
                "-",
                substring(samples$sample_date,5,6),
                "-",
                substring(samples$sample_date,7,8))
samples$Date <- as.Date(sample_date)
samples <- samples[,c(1,12,2:9,11)]

samplingDates <- unique(as.Date(sample_date))
samplingDates <- samplingDates[samplingDates > "2017-06-03" & samplingDates < "2017-08-24"] # Include only experiment
#samplingDates[samplingDates >= max(samplingDates)] <- "2017-08-23"

## Table 1 - All Chemistry 
# WC averaged over the entire experiment for each individual.
samplingDates_c <- paste0(substring(wc_c$sample_date,1,4),
                          "-",
                          substring(wc_c$sample_date,5,6),
                          "-",
                          substring(wc_c$sample_date,7,8))
wc_c$Date <- as.Date(samplingDates_c)
# Subset to important columns
wc_c_sub <- subset(wc_c,select=c("Date","tankID","temp","salinity","corr_CT","corr_AT","pH_out","pCO2_out","Ca_out","Ar_out"))

# Ugly loop that summarizes tank chemistry by sample time point
counter=1
for(j in unique(wc_c_sub$tankID)){
  wc_temp <- wc_c_sub[wc_c_sub$tankID == j,]
  assign <- NULL
  for(i in 1:length(samplingDates)){
    temp <- wc_temp$Date-samplingDates[i]# - wc_temp$Date
    temp2 <- ifelse(c(temp <= 0),1,0)
    assign <- cbind(assign,temp2)
  }
  colnames(assign) <- paste0("TP_",samplingDates-min(samplingDates)+1)
  for(i in 1:length(samplingDates)){
    if(i == 1 & counter == 1){
      wc_c_sum <- data.frame(Date=samplingDates[i],
                             tankID=j,
                             Temp_Complete=mean(wc_temp$temp[assign[,i]==1]),
                             Sal_Complete=mean(wc_temp$salinity[assign[,i]==1]),
                             CT_Complete=mean(wc_temp$corr_CT[assign[,i]==1]),
                             AT_Complete=mean(wc_temp$corr_AT[assign[,i]==1]),
                             pH_Complete=mean(wc_temp$pH_out[assign[,i]==1]),
                             pCO2_Complete=mean(wc_temp$pCO2_out[assign[,i]==1]),
                             CaSat_Complete=mean(wc_temp$Ca_out[assign[,i]==1]),
                             ArgSat_Complete=mean(wc_temp$Ar_out[assign[,i]==1]))
    }else{
      wc_c_sum <-rbind(wc_c_sum,
                       data.frame(Date=samplingDates[i],
                                  tankID=j,
                                  Temp_Complete=mean(wc_temp$temp[assign[,i]==1]),
                                  Sal_Complete=mean(wc_temp$salinity[assign[,i]==1]),
                                  CT_Complete=mean(wc_temp$corr_CT[assign[,i]==1]),
                                  AT_Complete=mean(wc_temp$corr_AT[assign[,i]==1]),
                                  pH_Complete=mean(wc_temp$pH_out[assign[,i]==1]),
                                  pCO2_Complete=mean(wc_temp$pCO2_out[assign[,i]==1]),
                                  CaSat_Complete=mean(wc_temp$Ca_out[assign[,i]==1]),
                                  ArgSat_Complete=mean(wc_temp$Ar_out[assign[,i]==1]))
                       )
    }
    counter = counter + 1
  }
}

## Table 2 - Two week Chemistry
# Weekly WC (temp, salinity, NBS pH) averaged based on the two previous weeks 
# prior to sample (used for delta pH)
wc_w$Date <- as.Date(wc_w$Date)
#write.csv(wc_w,"/home/downeyam/Github/AE17_Cvirginica_MolecularResponse/data/water_chem/AE17_weeklySummary.csv")

wc_w_sub <- subset(wc_w,select=c("Date","TankID","Temperature","Sal_Corr","pH_scaleFree","pH_NBS","pH_SW","pH_Total"))
# Remove pre exposure timepoints
## Ugle for loop that summarizes weekly water chem (temp, salinity, and pH) for each tank based on previous
 # 2 weeks

Interval <- 14 # number of days summarizing

for(j in unique(wc_w_sub$TankID)){
  wc_temp <- wc_w_sub[wc_w_sub$TankID == j,]
  assign <- NULL
  for(i in 1:length(samplingDates)){
    temp <- wc_temp$Date-samplingDates[i]
    temp2 <- ifelse(c(temp >= -Interval & temp <= 0),1,0)
    assign <- cbind(assign,temp2)
  }
  colnames(assign) <- paste0("TP_",samplingDates-min(samplingDates)+1)
  
  for(i in 1:length(samplingDates)){
    if(i == 1 & j == 1){
      wc_w_sum <- data.frame(Date=samplingDates[i],
                             tankID=j,
                             Temp_2W=mean(wc_temp$Temperature[assign[,i]==1]),
                             Sal_2W=mean(wc_temp$Sal_Corr[assign[,i]==1]),
                             pH_NBS_2W=mean(wc_temp$pH_NBS[assign[,i]==1]),
                             pH_SW_2W=mean(wc_temp$pH_SW[assign[,i]==1]),
                             pH_Total_2W=mean(wc_temp$pH_Total[assign[,i]==1]),
                             pH_scaleFree_2W=mean(wc_temp$pH_scaleFree[assign[,i]==1]))
    }else{
      wc_w_sum <-rbind(wc_w_sum,
                       data.frame(Date=samplingDates[i],
                                  tankID=j,
                                  Temp_2W=mean(wc_temp$Temperature[assign[,i]==1]),
                                  Sal_2W=mean(wc_temp$Sal_Corr[assign[,i]==1]),
                                  pH_NBS_2W=mean(wc_temp$pH_NBS[assign[,i]==1]),
                                  pH_SW_2W=mean(wc_temp$pH_SW[assign[,i]==1]),
                                  pH_Total_2W=mean(wc_temp$pH_Total[assign[,i]==1]),
                                  pH_scaleFree_2W=mean(wc_temp$pH_scaleFree[assign[,i]==1]))
      )
    }
  }
}

## Table 3 - Calcification Chemistry
# Summary of WC for the first 33 days of exposure. Used for calcification data comparison.
cal_date <- as.Date("2017-07-05") # the last WC time point before calcification measure

# Ugly loop that summarizes tank chemistry by sample time point
counter=1
for(j in unique(wc_c_sub$tankID)){
  wc_temp <- wc_c_sub[wc_c_sub$tankID == j,]
  #wc_temp <- wc_c_sub[wc_c_sub$tankID == 1,]
  temp <- wc_temp$Date-cal_date # - wc_temp$Date
  temp2 <- ifelse(c(temp <= 0),1,0)
  if(counter == 1 & j == 1){
        wc_cal_sum <- data.frame(Date=samplingDates[i],
                               tankID=j,
                               Temp_Cal=mean(wc_temp$temp[temp2==1]),
                               Sal_Cal=mean(wc_temp$salinity[temp2==1]),
                               CT_Cal=mean(wc_temp$corr_CT[temp2==1]),
                               AT_Cal=mean(wc_temp$corr_AT[temp2==1]),
                               pH_Cal=mean(wc_temp$pH_out[temp2==1]),
                               pCO2_Cal=mean(wc_temp$pCO2_out[temp2==1]),
                               CaSat_Cal=mean(wc_temp$Ca_out[temp2==1]),
                               ArgSat_Cal=mean(wc_temp$Ar_out[temp2==1]))
  }else{
        wc_cal_sum <-rbind(wc_cal_sum,
                         data.frame(Date=samplingDates[i],
                                    tankID=j,
                                    Temp_Cal=mean(wc_temp$temp[temp2==1]),
                                    Sal_Cal=mean(wc_temp$salinity[temp2==1]),
                                    CT_Cal=mean(wc_temp$corr_CT[temp2==1]),
                                    AT_Cal=mean(wc_temp$corr_AT[temp2==1]),
                                    pH_Cal=mean(wc_temp$pH_out[temp2==1]),
                                    pCO2_Cal=mean(wc_temp$pCO2_out[temp2==1]),
                                    CaSat_Cal=mean(wc_temp$Ca_out[temp2==1]),
                                    ArgSat_Cal=mean(wc_temp$Ar_out[temp2==1]))
        )
  }
  counter = counter + 1
}

### Merging WC data
wc_sum <- left_join(wc_c_sum,wc_w_sum,by=c("Date","tankID"))
write.csv(wc_sum,"data/water_chem/AE17_CompleteCarbChem.csv")

#### Merge phenotype / sample data ####
## EPF data
comb_epf <- left_join(samples,epf_sel,by="ID")
comb_epf <- inner_join(comb_epf,wc_sum,by=c("Date","tankID"))
comb_epf$delta_ph_NBS <-comb_epf$pHNBS-comb_epf$pH_NBS_2W
comb_epf$delta_ph_total <-comb_epf$pHTotal-comb_epf$pH_Total_2W
comb_epf_reduce <-  comb_epf[,c(1,12:31)]
write.csv(comb_epf,"data/Phenotype/AE17_EPFpHComplete.csv")

## Calcification Data
comb_cal <- inner_join(samples,cal_sel,by="ID")
comb_cal <- inner_join(comb_cal,comb_epf_reduce,by="ID")
comb_cal <- left_join(comb_cal,wc_cal_sum,by="tankID")
write.csv(comb_cal,"data/Phenotype/AE17_CalcificationComplete.csv")

