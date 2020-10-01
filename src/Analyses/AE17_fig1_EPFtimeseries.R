#### Script used to analyze the full EPF pH timeseries data ####
# Figures from this code : 1

# Note : pH evaluated on the total scale.

## packages
library(mgcv)
library(dplyr)
library(car)
library(lme4)
library(lmerTest)
library(factoextra)
library(multcomp)
library(multcompView)
library("RColorBrewer")
pal <- brewer.pal(n = 12, name = 'Paired')
yellow <- brewer.pal(n=9,name = 'YlOrRd')[5]
col_perm <- c(pal[1:2],pal[5:6],yellow)
#col_perm <- c(pal[1:2],pal[5:6],pal[12])
# Located in src Analysis/Phenotype folder, will need to set full working directory or setwd()
setwd("/home/downeyam/Github/AE17_Cvirginica_MolecularResponse")
source("src/Accessory/basicR_functions.R")

#### Sample Data ####
epf_exp <- read.csv("data/Phenotype/CompletePhenotype.csv",stringsAsFactors = FALSE)
epf_exp$EPF_pH <- epf_exp$EPF_pH_Total
epf_exp$pCO2_fac <- as.factor(epf_exp$pCO2)
epf_exp$Timepoint_fac <- as.factor(epf_exp$Timepoint)
epf_exp$EPF_envAdj <- epf_exp$EPF_pH-epf_exp$pH_Total_2W

#### Water Chemistry Data ####
wc <- read.delim("data/water_chem/AE17_WaterChemistry_weekly.csv",sep=",",
                 stringsAsFactors = FALSE)


#### Analysis ####
# Random effects model using lmer

#### Measured pH vs. treatment and time ####
## Full Model 
epfAllTP_full <- lmer(EPF_pH~pCO2_fac*Timepoint_fac + (1|PopOrigin) + (1|TankID),data=epf_exp) 

## Reduced (final) model 
epfAllTP_red <- lmer(EPF_pH~pCO2_fac*Timepoint_fac + (1|TankID),data=epf_exp)
# LRT to check if reduced model is diff. than full model
anova(epfAllTP_full,epfAllTP_red) # It isnt, so we use the simpler reduced model

## Check assumptions
plot(epfAllTP_red)

anova(epfAllTP_red)
#Type III Analysis of Variance Table with Satterthwaite's method
#                        Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)   
#pCO2_fac               0.80996 0.40498     2 15.605  6.4073 0.009302 **
#Timepoint_fac          0.17505 0.03501     5 85.947  0.5539 0.734949   
#pCO2_fac:Timepoint_fac 1.72395 0.17240    10 85.931  2.7275 0.005833 **

# Planned Comparisons  
# Create 'k' object with planned comparisons.
# These will include comparisons within each timepoint for
# for either OA treatment (900 or 2800) vs the control.
group <- paste0(epf_exp$pCO2_fac,":",epf_exp$Timepoint_fac)
mod_matrix <- model.matrix(epfAllTP_red)
agg_mod_matrix <- aggregate(mod_matrix~group,FUN=mean)
rownames(agg_mod_matrix) <- agg_mod_matrix$group
agg_mod_matrix <- agg_mod_matrix[,-1]
lc2 <- as.matrix(agg_mod_matrix)
k <- rbind("400v900_1"=lc2["900:1",]    -lc2["400:1",],
           "400v2800_1"=lc2["2800:1",]  -lc2["400:1",],
           "400v900_2"=lc2["900:2",]    -lc2["400:2",],
           "400v2800_2"=lc2["2800:2",]  -lc2["400:2",],
           "400v900_9"=lc2["900:9",]    -lc2["400:9",],
           "400v2800_9"=lc2["2800:9",]  -lc2["400:9",],
           "400v900_22"=lc2["900:22",]  -lc2["400:22",],
           "400v2800_22"=lc2["2800:22",]-lc2["400:22",],
           "400v900_50"=lc2["900:50",]  -lc2["400:50",],
           "400v2800_50"=lc2["2800:50",]-lc2["400:50",],
           "400v900_80"=lc2["900:80",]  -lc2["400:80",],
           "400v2800_80"=lc2["2800:80",]-lc2["400:80",]
)
(epfAllTP_posthoc_Correction <- summary(glht(epfAllTP_red,linfct=k),adjusted(type = "fdr")))

#### Relative pH vs. treatment and time ####
## Full Model
epfenvAllTP_full <- lmer(EPF_envAdj~pCO2_fac*Timepoint_fac + (1|PopOrigin) + (1|TankID),data=epf_exp) 
## Reduced (final) Model
epfenvAllTP_red <- lmer(EPF_envAdj~pCO2_fac*Timepoint_fac + (1|TankID),data=epf_exp)
# LRT to check if reduced model is diff. than full model
anova(epfenvAllTP_red,epfenvAllTP_full) # It isnt, so we use the simpler reduced model
# ANOVA 
anova(epfenvAllTP_red)
#Type III Analysis of Variance Table with Satterthwaite's method
#                        Sum Sq Mean Sq NumDF DenDF F value  Pr(>F)  
#pCO2_fac               0.78091 0.39045     2    12  6.4769 0.01237 *
#Timepoint_fac          0.55850 0.11170     5    60  1.8529 0.11621  
#pCO2_fac:Timepoint_fac 1.19071 0.11907    10    60  1.9752 0.05220 .

# Planned comparisons (using same comparisons when looking at measured pH as response variable)
(epfAllTP_posthoc_Correction <- summary(glht(epfenvAllTP_red,linfct=k),adjusted(type = "fdr")))

## T test to look at delta EPF relative to seawater pH
tTest_out <- NULL
#x<-t.test(epf_exp$EPF_envAdj[epf_exp$Timepoint == unique(epf_exp$Timepoint)[1] & epf_exp$pCO2 == unique(epf_exp$pCO2)[1]],mu = 0)
x<-NULL
for(i in unique(epf_exp$Timepoint)){
   for(j in unique(epf_exp$pCO2)){
      #print(j)
      x<-c(x,t.test(epf_exp$EPF_envAdj[epf_exp$Timepoint == i & epf_exp$pCO2 == j],
                    mu = 0)$p.value)
   }
}

y <- expand.grid(unique(epf_exp$pCO2),unique(epf_exp$Timepoint))
tTest_df <- data.frame(Timepoint=y[,2],pCO2=y[,1],pValue=x)
tTest_df$pValue_corr <- p.adjust(tTest_df$pValue,method = "BH")
tTest_df$Sig <- ifelse(tTest_df$pValue_corr < 0.05,"TRUE","FALSE")
tTest_df <- tTest_df[order(tTest_df$Timepoint),]
tTest_df <- tTest_df[order(tTest_df$pCO2),]
x <- NULL
for(j in 1:nrow(tTest_df)){
   i <- tTest_df$pValue_corr[j]
   if(i >= 0.1){x <- c(x,"")}
   if(i < 0.1 & i >= 0.05){x <- c(x,"+")}
   if(i < 0.05 & i >= 0.01){x <- c(x,"*")}
   if(i < 0.01 & i >= 0.001){x <- c(x,"**")}
   if(i < 0.001){x <- c(x,"***")}
}
tTest_df$sig_signal <- x

# Quick visualization of delta pH to make sure tTest values make sense
ggplot(epf_exp,aes(y=EPF_envAdj,x=Timepoint_fac,colour=pCO2_fac)) + 
   geom_boxplot() + geom_hline(yintercept = 0) +
   geom_jitter(width = 0.2) +
   facet_grid(.~pCO2) 

#### Figure #####################   ##################################################
# This figure uses base plot() to create a two-panel figure with the EPF timeseries
## Summarize the pH of the environment (featured as lines on figure)
c_mean <- mean(wc$pH_Total[wc$PCO2 == 550])
oa_900_mean <- mean(wc$pH_Total[wc$PCO2 == 1000])
oa_2800_mean <- mean(wc$pH_Total[wc$PCO2 == 2800])
#oa_2800_mean <- mean(epf_exp$pH_Complete[epf_exp$pCO2_fac == 2800])
# Code for setting the plotting space for the two panels (plus middle section for shared legend)
m <- matrix(c(1,2,3),nrow = 3,ncol = 1,byrow = TRUE)
layout(mat = m,heights = c(0.5,0.05,0.5),widths = c(.1))
par(mar = c(5,5,3,2))

#### Panel A - measured EPF pH ####
### Summarize measured EPF pH for plotting
## Take the means for each timepoint and treatment level
treatSeq_means <- aggregate(EPF_pH~Timepoint+pCO2_fac,epf_exp,FUN=mean)
treatSeq_means$pCO2_fac <- as.character(treatSeq_means$pCO2_fac)
treatSeq_means$pCO2_fac[treatSeq_means$pCO2_fac == "400"]  <- "Control"
treatSeq_means$pCO2_fac[treatSeq_means$pCO2_fac == "900"]  <- "Mod. OA"
treatSeq_means$pCO2_fac[treatSeq_means$pCO2_fac == "2800"]  <- "High OA"
treatSeq_means$pCO2_fac <- as.factor(treatSeq_means$pCO2_fac )
treatSeq_means$pCO2_fac <- factor(treatSeq_means$pCO2_fac, levels = c("Control", "Mod. OA", "High OA"))
## Take the standard error (SE) for each timepoint and treatment level
treatSeq_SE <- aggregate(EPF_pH~Timepoint+pCO2_fac,epf_exp,FUN=se)
treatSeq_means$error <- treatSeq_SE$EPF_pH
library(cowplot)

# EPF pH plot
pA <- ggplot(treatSeq_means,aes(x=Timepoint,y=EPF_pH,group=pCO2_fac,shape=pCO2_fac,colour=pCO2_fac)) +
   geom_hline(yintercept=c_mean,colour=col_perm[2],linetype=17,size=.8,show.legend =TRUE) + # control 
   geom_hline(yintercept=oa_900_mean,colour=col_perm[5],linetype=16,size=.8) +# 900
   geom_hline(yintercept=oa_2800_mean,colour=col_perm[4],linetype=15,size=.8) +# 2800
   geom_line(position=position_dodge(width=0.15)) +
   geom_point(size=5,position=position_dodge(width=0.15)) +
   geom_errorbar(aes(ymin=EPF_pH-error,ymax=EPF_pH+error),width=0.1,position=position_dodge(width=0.15)) +
   scale_y_continuous(breaks=c(6.50,7.00,7.50,8.00),limits=c(6.45,8.15),labels=c(" 6.50"," 7.00"," 7.50"," 8.00")) +
   scale_x_log10(breaks = c(1,2,9,22,50,80)) +
   scale_shape_manual(values=c(16,15,17))+
   scale_colour_manual(values=c(col_perm[2],col_perm[5],col_perm[4])) +
   labs(x="Time (Days)",y=expression(pH[EPF]~(Total)),fill="") +
   theme_cowplot() +
   theme(panel.border = element_blank(),
      legend.direction = "horizontal",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      legend.spacing.y = unit(10, "mm"),
      legend.position = c(0.02,.08),
      legend.background = element_rect(linetype = 1, 
                                       size = 0.5, 
                                       colour = 1),
      legend.margin = margin(6, 6, 6, 6))
pA

# Delta pH plot
treatSeq_means <- aggregate(EPF_envAdj~Timepoint+pCO2_fac,epf_exp,FUN=mean)
treatSeq_means$pCO2_fac <- as.character(treatSeq_means$pCO2_fac)
treatSeq_means$pCO2_fac[treatSeq_means$pCO2_fac == "400"]  <- "Control"
treatSeq_means$pCO2_fac[treatSeq_means$pCO2_fac == "900"]  <- "Mod. OA"
treatSeq_means$pCO2_fac[treatSeq_means$pCO2_fac == "2800"]  <- "High OA"
treatSeq_means$pCO2_fac <- as.factor(treatSeq_means$pCO2_fac )
treatSeq_means$pCO2_fac <- factor(treatSeq_means$pCO2_fac, levels = c("Control", "Mod. OA", "High OA"))
treatSeq_SE <- aggregate(EPF_envAdj~Timepoint+pCO2_fac,epf_exp,FUN=se)
treatSeq_means$error <- treatSeq_SE$EPF_envAdj
pB <- ggplot(treatSeq_means,aes(x=Timepoint,y=EPF_envAdj,group=pCO2_fac,shape=pCO2_fac,colour=pCO2_fac)) +
   #geom_hline(yintercept=0,linetype=3,size=0.8) +
   geom_line(position=position_dodge(width=0.15)) +
   geom_point(size=5,position=position_dodge(width=0.15)) +
   geom_errorbar(aes(ymin=EPF_envAdj-error,ymax=EPF_envAdj+error),width=0.1,position=position_dodge(width=0.15)) +
   scale_y_continuous(limits=c(-1.1,0.3),breaks=c(-0.9,-0.60,-0.30,0,0.3),labels=c("-0.90","-0.60","-0.30","0.00","0.30")) +
   scale_x_log10(breaks = c(1,2,9,22,50,80)) +
   scale_shape_manual(values=c(16,15,17))+
   scale_colour_manual(values=c(col_perm[2],col_perm[5],col_perm[4])) +
   labs(x="Time (Days)",y=expression(paste(Delta," pH (Total)")),fill="") +
   theme_cowplot() +
   theme(panel.border = element_blank(),
         legend.direction = "horizontal",
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.title = element_blank(),
         #legend.spacing.y = unit(0, "mm"),
         legend.position = c(0.02,.08),
         legend.background = element_rect(linetype = 1, 
                                          size = 0.5, 
                                          colour = 1),
         legend.margin = margin(6, 6, 6, 6)) 

plot_grid("",pA,"",pB,labels = c("A","","B",""),nrow=4,rel_heights = c(1,5,2,5))
