## packages
library(mgcv)
library(dplyr)
library(car)
library(lme4)
library(lmerTest)
library(factoextra)
library(multcomp)
library(multcompView)
library(cowplot)
library(grid)
library(gridExtra)
library("RColorBrewer")
pal <- brewer.pal(n = 12, name = 'Paired')
col_perm <- c(pal[1:2],pal[5:6],pal[12])
# Located in src Analysis/Phenotype folder, will need to set full working directory or setwd()
setwd("/home/downeyam/Github/AE17_Cvirginica_MolecularResponse/src/Analyses/Phenotype/")
source("basicR_functions.R")

#### Data ####
setwd("/home/downeyam/Github/AE17_Cvirginica_MolecularResponse/data/Phenotype/")
pheno <- readRDS("AE17_summaryPhenotype_exposure.RData") # just exposure timepoints
pheno2 <- readRDS("AE17_summaryPhenotype_alltimepoints.RData") # acclimation and exposure timepoints
# Remove 81 which we know has wonky EPF values (for consistency these are also ultimately removed from the calcification estimates as well)
pheno_red <- pheno[pheno$timepoint != 81,]

epf_exp <- pheno_red[!is.na(pheno_red$EPF_pH),]
epf_exp <-epf_exp[as.numeric(epf_exp$timepoint) > 40,]

## Full Model
epfAllTP_full <- lmer(EPF_pH~pCO2_fac*timepoint_fac + (1|PopOrigin) + (1|shelf/tank),data=epf_exp) 
## Reduced (final) model 
epfAllTP_red <- aov(EPF_pH~pCO2_fac,data=epf_exp)
#epfAllTP_red <- aov(EPF_pH~pCO2_fac*timepoint_fac,data=epf_exp)
# ANOVA
summary(epfAllTP_red)
# Checking reduced model is an improvement
anova(epfAllTP_red,epfAllTP_full)

TukeyHSD(epfAllTP_red)

#### Figure ####

c_mean <- mean(epf_exp$Tank_pH[epf_exp$pCO2_fac==400])
oa_900_mean <- mean(epf_exp$Tank_pH[epf_exp$pCO2_fac==900])
oa_2800_mean <- mean(epf_exp$Tank_pH[epf_exp$pCO2_fac==2800])


## Figures

## Measured EPF ph plot

# Bar plot
treatSeq_means <- aggregate(EPF_pH~pCO2_fac,epf_exp,FUN=mean)
treat_means_ctrl <- treatSeq_means[treatSeq_means$pCO2_fac == "400",]
treat_means_oa_900 <- treatSeq_means[treatSeq_means$pCO2_fac == "900",]
treat_means_oa_2800 <- treatSeq_means[treatSeq_means$pCO2_fac == "2800",]

treatSeq_SE <- aggregate(EPF_pH~pCO2_fac,epf_exp,FUN=se)
treat_SE_ctrl <- treatSeq_SE[treatSeq_SE$pCO2_fac == "400",]
treat_SE_oa_900 <- treatSeq_SE[treatSeq_SE$pCO2_fac == "900",]
treat_SE_oa_2800 <- treatSeq_SE[treatSeq_SE$pCO2_fac == "2800",]

c_mean <- mean(epf_exp$Tank_pH[epf_exp$pCO2_fac==400])
oa_900_mean <- mean(epf_exp$Tank_pH[epf_exp$pCO2_fac==900])
oa_2800_mean <- mean(epf_exp$Tank_pH[epf_exp$pCO2_fac==2800])

bar_measuredpH <- ggplot(treatSeq_means,aes(x=pCO2_fac,
                                            y=EPF_pH,
                                            fill=pCO2_fac)) + 
  # Mean environment lines underneath other lines   
  geom_hline(yintercept=c_mean,colour=col_perm[2],size=.8,linetype=2) + 
  geom_hline(yintercept=oa_900_mean,colour=col_perm[5],size=.8,linetype=2) + 
  geom_hline(yintercept=oa_2800_mean,colour=col_perm[4],size=.8,linetype=2) + 
  geom_bar(stat="identity",width=0.5) +
  geom_errorbar(aes(ymin=treatSeq_means$EPF_pH-treatSeq_SE$EPF_pH,
                    ymax=treatSeq_means$EPF_pH+treatSeq_SE$EPF_pH), 
                width=0.1,size=0.8) +
  coord_cartesian(ylim=c(6,8)) + 
  theme_cowplot() + 
  scale_fill_manual(values=c(col_perm[2],col_perm[5],col_perm[4]))+
  labs(x="",y="EPF pH (NBS)",fill="") +
  theme(legend.position = "NULL")
bar_measuredpH

# Box plot
box_measuredpH <- ggplot(epf_exp,aes(x=pCO2_fac,y=EPF_pH,fill=pCO2_fac)) + 
  geom_hline(yintercept=c_mean,colour=col_perm[2],size=.8,linetype=2) + 
  geom_hline(yintercept=oa_900_mean,colour=col_perm[5],size=.8,linetype=2) + 
  geom_hline(yintercept=oa_2800_mean,colour=col_perm[4],size=.8,linetype=2) +
  geom_boxplot() +
  theme_cowplot() + 
  scale_fill_manual(values=c(col_perm[2],col_perm[5],col_perm[4]))+
  labs(x="",y="EPF pH (NBS)",fill="") +
  theme(legend.position = "NULL")
box_measuredpH 

## Relative EPF pH

#Bar plot 
# Bar plot
treatRSeq_means <- aggregate(EPF_envAdj~pCO2_fac,epf_exp,FUN=mean)
treatR_means_ctrl <- treatRSeq_means[treatRSeq_means$pCO2_fac == "400",]
treatR_means_oa_900 <- treatRSeq_means[treatRSeq_means$pCO2_fac == "900",]
treatR_means_oa_2800 <- treatRSeq_means[treatRSeq_means$pCO2_fac == "2800",]

treatRSeq_SE <- aggregate(EPF_envAdj~pCO2_fac,epf_exp,FUN=se)
treatR_SE_ctrl <- treatRSeq_SE[treatRSeq_SE$pCO2_fac == "400",]
treatR_SE_oa_900 <- treatRSeq_SE[treatRSeq_SE$pCO2_fac == "900",]
treatR_SE_oa_2800 <- treatRSeq_SE[treatRSeq_SE$pCO2_fac == "2800",]

bar_relativepH <- ggplot(treatRSeq_means,aes(x=pCO2_fac,
                                            y=EPF_envAdj,
                                            fill=pCO2_fac)) + 
  # Mean environment lines underneath other lines   
  geom_hline(yintercept=0,linetype=2,size=0.8) +
  geom_bar(stat="identity",width=0.5) +
  geom_errorbar(aes(ymin=treatRSeq_means$EPF_envAdj-treatRSeq_SE$EPF_envAdj,
                    ymax=treatRSeq_means$EPF_envAdj+treatRSeq_SE$EPF_envAdj), 
                width=0.1,size=0.8) +
  coord_cartesian(ylim=c(-0.6,0.6)) + 
  theme_cowplot() + 
  scale_fill_manual(values=c(col_perm[2],col_perm[5],col_perm[4]))+
  labs(x="",y="EPF pH - Environment pH(NBS)",fill="") +
  theme(legend.position= "NULL")
bar_relativepH

# Box plot
box_relativepH <- ggplot(epf_exp,aes(x=pCO2_fac,y=EPF_envAdj,fill=pCO2_fac)) + 
  geom_hline(yintercept=0,linetype=2,size=0.8) +
  geom_boxplot() +
  theme_cowplot() + 
  scale_fill_manual(values=c(col_perm[2],col_perm[5],col_perm[4]))+
  labs(x="",y="EPF pH - Environment pH(NBS)",fill="") +
  theme(legend.position = "NULL")
box_relativepH

x.grob1 <- textGrob("pCO2 (uatm)", 
                    gp=gpar(col="black", fontsize=15))

boxes <- plot_grid(box_measuredpH,box_relativepH,labels=c("A","B"))
grid.arrange(arrangeGrob(boxes,bottom = x.grob1))

bars <- plot_grid(bar_measuredpH,bar_relativepH,labels=c("A","B"))
grid.arrange(arrangeGrob(bars,bottom = x.grob1))

