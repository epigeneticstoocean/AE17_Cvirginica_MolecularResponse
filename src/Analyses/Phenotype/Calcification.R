#### AE17 Calcification Analysis ####

## packages
knitr::opts_chunk$set(echo = TRUE)
library(dplyr)
library(ggplot2)
library(lmerTest)
library(car)
library(RColorBrewer)
pal <- brewer.pal(n = 12, name = 'Paired')
col_perm <- c(pal[1:2],pal[5:6],pal[12])
# Located in src Analysis/Phenotype folder, will need to set full working directory or setwd()
source("basicR_functions.R")

#### Data ####
setwd("/home/downeyam/Github/AE17_Cvirginica_MolecularResponse/data/Phenotype/")
# Final version of bouyant weights
cal<-read.csv("AE17_BouyantWeight_Final.csv")
# Final version of sample data
samp<-read.csv("AE17_CollectionInfo.csv")
# Final EPF fluid data 
pheno <- readRDS("AE17_summaryPhenotype_exposure.RData")
pheno_red <- pheno[pheno$timepoint != 81,]
pheno_red <- subset(pheno_red,select=c(1,30,37))

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

#### Analysis #####

#### Calcification vs Environment ####

## Full model
cal_full <- lmer(calcification~pCO2_calc + (1|PopOrigin) + (1|shelf/tank),data=cal_red2)
## Reduced model final)
cal_analysis_fixed <- lm(calcification~pCO2_calc,data=cal_red2)


#### Calcification vs EPF pH ####
cal_epf <- inner_join(cal_red2,pheno_red)
# Full Model
calVspH_full <- lmer(calcification ~ EPF_pH+sample_date + (1|PopOrigin) + (1|shelf/tank),data=cal_epf)
# Reducded Model (final)
calVspH_simple <- lm(calcification ~ EPF_pH,data=cal_epf)

#### Figure ####
## Plot 1 
#Mean calcification rate 
mean_cal <- aggregate(calcification*100~pCO2,cal_red2,FUN=mean)
se_cal <- aggregate(calcification*100~pCO2,cal_red2,FUN=se)
se_pco2 <- aggregate(pCO2_calc~pCO2,cal_red2,FUN=sd)
mean_cal <- data.frame(pCO2=mean_cal$pCO2,Rel_Change=mean_cal$calcification,
                       ymin=mean_cal$calcification-c(se_cal$calcification*1.96),
                       ymax=mean_cal$calcification+c(se_cal$calcification*1.96),
                       xmin=mean_cal$pCO2-c(se_pco2$pCO2_calc*1.96),
                       xmax=mean_cal$pCO2+c(se_pco2$pCO2_calc*1.96))

# ggplot version
out <- summary(lm(calcification*100~pCO2,data=cal_red2))

p <- ggplot(mean_cal,aes(x=pCO2,y=Rel_Change,shape=as.factor(pCO2),colour=as.factor(pCO2))) +
  geom_abline(slope = out$coefficients[2,1],intercept = out$coefficients[1,1]) +
  geom_hline(aes(yintercept=0),linetype="dotted") +
  geom_point(aes(size=1.5)) +
  ylim(-0.05,0.05) + 
  xlim(300,3200) +
  scale_shape_manual(values=c(16,15,17))+
  scale_colour_manual(values=c(col_perm[2],col_perm[5],col_perm[4]))+
  geom_errorbarh(aes(xmin=xmin, xmax=xmax)) + 
  geom_errorbar(aes(ymin=ymin,ymax=ymax),width=75) 

t <- p + theme_bw(base_size = 16) + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p1<-t + ylab("Calcification (% change per day)") + xlab("pCO2 (uatm)") +
  coord_cartesian(clip = 'off') +
  geom_text(
    x = 100,
    y = 0.06,
    inherit.aes = FALSE,
    label = "A",
    check_overlap = FALSE,
    hjust = 1,
    size = 10
  ) +
  theme(legend.position="none",plot.margin = unit(c(5, 1, 1, 1), "lines")) 
p1

## Plot 2
cal_epf$pCO2_name <- "NA"
cal_epf$pCO2_name[cal_epf$pCO2 == unique(cal_epf$pCO2)[1]] <-  "Ambient"
cal_epf$pCO2_name[cal_epf$pCO2 == unique(cal_epf$pCO2)[2]] <-  "OA 1000" 
cal_epf$pCO2_name[cal_epf$pCO2 == unique(cal_epf$pCO2)[3]] <-  "OA 2800"
cal_epf$pCO2_name <- factor(cal_epf$pCO2_name,levels = c("Ambient", "OA 1000", "OA 2800"))

p <- ggplot(cal_epf,aes(x=EPF_pH,y=calcification*100,colour=pCO2_name)) +
  geom_point(size=4,aes(shape=pCO2_name)) + 
  ylim(-0.05,0.05) +
  scale_color_manual(values=c(col_perm[2],col_perm[5],col_perm[4])) +
  scale_shape_manual(values=c(16,15,17)) +
  geom_abline(slope = calVspH_simple$coefficients[2],
              intercept = calVspH_simple$coefficients[1]) 

t <- p + theme_bw(base_size = 16) + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.15,.85),
        axis.line = element_line(colour = "black")) 

p2<-t + ylab("Calcification (% change per day)") + xlab("EPF pH (NBS)") + labs(colour="pCO2 (uatm)",shape="pCO2 (uatm)") +
  guides(shape = guide_legend(override.aes = list(size = 4))) +
  coord_cartesian(clip = 'off') +
  geom_text(
    x = 6.7,
    y = 0.06,
    inherit.aes = FALSE,
    label = "B",
    check_overlap = FALSE,
    hjust = 1,
    size = 10
  ) +
  theme(plot.margin = unit(c(5, 1, 1, 1), "lines"))

multiplot(p1,p2,cols=2)
