#### Script used to analyze the full EPF pH timeseries data ####
# Figures from this code : 1

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
col_perm <- c(pal[1:2],pal[5:6],pal[12])
# Located in src Analysis/Phenotype folder, will need to set full working directory or setwd()
setwd("/home/downeyam/Github/AE17_Cvirginica_MolecularResponse")
source("src/Accessory/basicR_functions.R")

#### Sample Data ####
epf_exp <- read.csv("data/Phenotype/CompletePhenotype_final2020.csv",stringsAsFactors = FALSE)
epf_exp$EPF_pH <- epf_exp$EPF_Total_pH
epf_exp$pCO2_fac <- as.factor(epf_exp$pCO2)
epf_exp$Timepoint_fac <- as.factor(epf_exp$Timepoint)
epf_exp$EPF_envAdj <- epf_exp$EPF_pH-epf_exp$pH_Total_2W
epf_exp$EPF_envAdjAlt <- epf_exp$EPF_pH-epf_exp$pH_scaleFree_2W

# Checking difference between NBS and scale free pH
plot(epf_exp$EPF_envAdj~epf_exp$EPF_envAdjAlt)
abline(a=0,b=1,col="red")

#### Water Chemistry Data ####
wc <- read.delim("data/water_chem/AE17_weeklyExposure_final2020.csv",sep=",",
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
#                        Sum Sq  Mean Sq NumDF  DenDF F value  Pr(>F)  
#pCO2_fac               0.58311 0.291553     2 14.854  4.9360 0.02272 *
#Timepoint_fac          0.27384 0.054769     5 73.993  0.9272 0.46837  
#pCO2_fac:Timepoint_fac 1.40569 0.140569    10 73.985  2.3798 0.01666 *

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
x<-t.test(epf_exp$EPF_envAdj[epf_exp$Timepoint == unique(epf_exp$Timepoint)[1] & epf_exp$pCO2 == unique(epf_exp$pCO2)[1]],mu = 0)
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
treatSeq_means$Timepoint[treatSeq_means$Timepoint == 2] <- 4
treat_means_ctrl <- treatSeq_means[treatSeq_means$pCO2_fac == "400",]
treat_means_oa_900 <- treatSeq_means[treatSeq_means$pCO2_fac == "900",]
treat_means_oa_2800 <- treatSeq_means[treatSeq_means$pCO2_fac == "2800",]
# Stagger the timepoints for 900 and 2800 treatments so they can be seen better on plot
treat_means_oa_900$timepoint <- treat_means_oa_900$Timepoint + 1
treat_means_oa_2800$timepoint <- treat_means_oa_2800$Timepoint + 2
# Manually alter timepoint 2 placement along x-axis for visual clarity
treat_means_oa_900$Timepoint <- treat_means_oa_900$Timepoint + 1
treat_means_oa_2800$Timepoint <- treat_means_oa_2800$Timepoint +2
## Take the standard error (SE) for each timepoint and treatment level
treatSeq_SE <- aggregate(EPF_pH~Timepoint+pCO2_fac,epf_exp,FUN=se)
treat_SE_ctrl <- treatSeq_SE[treatSeq_SE$pCO2_fac == "400",]$EPF_pH
treat_SE_oa_900 <- treatSeq_SE[treatSeq_SE$pCO2_fac == "900",]$EPF_pH
treat_SE_oa_2800 <- treatSeq_SE[treatSeq_SE$pCO2_fac == "2800",]$EPF_pH

# Base plot with 400
bp <- plot(treat_means_ctrl$EPF_pH~treat_means_ctrl$Timepoint,
           ylab=expression(pH[EPF]~(NBS)),xlab="Time (Days)",
           col=col_perm[2],pch=16,cex=2,cex.axis = 1.5,cex.lab=1.5,
           ylim = c(7.0,8.25),xlim=c(0,82),bty="n")

# Mean environment lines underneath other lines    
abline(h=c_mean,col=col_perm[2],lty=2) # control
abline(h=oa_900_mean,col=col_perm[5],lty=2) # 900
abline(h=oa_2800_mean,col=col_perm[4],lty=2) # 2800

# 400 EPF lines
lines(treat_means_ctrl$EPF_pH~treat_means_ctrl$Timepoint,
      col="lightblue4")
arrows(x0 = treat_means_ctrl$Timepoint,
       x1 = treat_means_ctrl$Timepoint,
       y0 = c(treat_means_ctrl$EPF_pH - treat_SE_ctrl),
       y1 = treat_means_ctrl$EPF_pH + treat_SE_ctrl,
       angle = 90, len = 0.05,
       code = 3, xpd = NA, lwd = 2,
       col=col_perm[2])
# 900 EPF lines
points(treat_means_oa_900$EPF_pH~treat_means_oa_900$Timepoint,
       col=col_perm[5],pch=15,cex=2)
lines(treat_means_oa_900$EPF_pH~treat_means_oa_900$Timepoint,
      col=col_perm[5])
arrows(x0 = treat_means_oa_900$Timepoint,
       x1 = treat_means_oa_900$Timepoint,
       y0 = c(treat_means_oa_900$EPF_pH - treat_SE_ctrl),
       y1 = treat_means_oa_900$EPF_pH + treat_SE_ctrl,
       angle = 90, len = 0.05,
       code = 3, xpd = NA, lwd = 2,
       col=col_perm[5])
# 2800 EPF lines 
points(treat_means_oa_2800$EPF_pH~treat_means_oa_2800$Timepoint,
       col=col_perm[4],pch=17,cex=2)
lines(treat_means_oa_2800$EPF_pH~treat_means_oa_2800$Timepoint,
      col=col_perm[4])
arrows(x0 = treat_means_oa_2800$Timepoint,
       x1 = treat_means_oa_2800$Timepoint,
       y0 = c(treat_means_oa_2800$EPF_pH - treat_SE_oa_2800),
       y1 = treat_means_oa_2800$EPF_pH + treat_SE_oa_2800,
       angle = 90, len = 0.05, code = 3, xpd = NA, lwd = 2,
       col=col_perm[4])

# Based on planned comparisons, significance stars
text(x=5+1,y=8.25,label="*",cex = 2.2,col=col_perm[4], xpd = NA)
text(x=51+1,y=8.25,label="*",cex = 2.2,col=col_perm[4], xpd = NA)
text(x=80+1,y=8.25,label="*",cex = 2.2,col=col_perm[4], xpd = NA)

# Letter labels
text(x=-8,y=8.35,label="A",cex = 2.2, xpd = NA)

#### Legend #####
par(mar = c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend("bottom",horiz = TRUE,
       legend = c("Control","Moderate OA","High OA"),
       col = c(col_perm[2],col_perm[5],col_perm[4]),
       pch=c(16,15,17),
       cex = 1.5,
       lwd = 2,
       bty = "n")
#### Panel B - diff EPF #### 
par(mar = c(5,5,3,2))
## Take the means for each timepoint and treatment level
treatSeq_means <- aggregate(EPF_envAdj~Timepoint+pCO2_fac,epf_exp,FUN=mean)
treatSeq_means$Timepoint[treatSeq_means$Timepoint == 2] <- 4 
treat_means_ctrl <- treatSeq_means[treatSeq_means$pCO2_fac == "400",]
treat_means_oa_900 <- treatSeq_means[treatSeq_means$pCO2_fac == "900",]
treat_means_oa_2800 <- treatSeq_means[treatSeq_means$pCO2_fac == "2800",]
# Stagger the timepoints for 900 and 2800 treatments so they can be seen better on plot
# Manually alter timepoint 2 placement along x-axis for visual clarity
treat_means_oa_900$Timepoint <- treat_means_oa_900$Timepoint+ 1
treat_means_oa_2800$Timepoint <- treat_means_oa_2800$Timepoint+ 2
## Take the standard error (se) for each timepoint and treatment level
treatSeq_SE <- aggregate(EPF_envAdj~Timepoint+pCO2_fac,epf_exp,FUN=se)
treat_SE_ctrl <- treatSeq_SE[treatSeq_SE$pCO2_fac == "400",]$EPF_envAdj
treat_SE_oa_900 <- treatSeq_SE[treatSeq_SE$pCO2_fac == "900",]$EPF_envAdj
treat_SE_oa_2800 <- treatSeq_SE[treatSeq_SE$pCO2_fac == "2800",]$EPF_envAdj

# Base plot with 400
bp <- plot(treat_means_ctrl$EPF_envAdj~treat_means_ctrl$Timepoint,
           ylab=expression(paste(Delta," pH (NBS)")),
           xlab="Time (Days)",pch=16,
           col=col_perm[2],cex=2,cex.axis = 1.5,cex.lab=1.5,
           ylim = c(-1.2,1.1),xlim=c(0,82),bty="n")
# Mean environment lines underneath other lines    
abline(h=0,col="black",lty=2) # control
# 400 EPF lines
lines(treat_means_ctrl$EPF_envAdj~treat_means_ctrl$Timepoint,
      col=col_perm[2])
arrows(x0 = treat_means_ctrl$Timepoint,
       x1 = treat_means_ctrl$Timepoint,
       y0 = c(treat_means_ctrl$EPF_envAdj - treat_SE_ctrl),
       y1 = treat_means_ctrl$EPF_envAdj + treat_SE_ctrl,
       angle = 90, len = 0.05,
       code = 3, xpd = NA, lwd = 2,
       col=col_perm[2])
# 900 EPF lines
points(treat_means_oa_900$EPF_envAdj~treat_means_oa_900$Timepoint,
       col=col_perm[5],pch=15,cex=2)
lines(treat_means_oa_900$EPF_envAdj~treat_means_oa_900$Timepoint,
      col=col_perm[5])
arrows(x0 = treat_means_oa_900$Timepoint,
       x1 = treat_means_oa_900$Timepoint,
       y0 = c(treat_means_oa_900$EPF_envAdj - treat_SE_oa_900),
       y1 = treat_means_oa_900$EPF_envAdj + treat_SE_oa_900,
       angle = 90, len = 0.05,
       code = 3, xpd = NA, lwd = 2,
       col=col_perm[5])
# 2800 EPF lines 
points(treat_means_oa_2800$EPF_envAdj~treat_means_oa_2800$Timepoint,
       col=col_perm[4],pch=17,cex=2)
lines(treat_means_oa_2800$EPF_envAdj~treat_means_oa_2800$Timepoint,
      col=col_perm[4])
arrows(x0 = treat_means_oa_2800$Timepoint,
       x1 = treat_means_oa_2800$Timepoint,
       y0 = c(treat_means_oa_2800$EPF_envAdj - treat_SE_oa_2800),
       y1 = treat_means_oa_2800$EPF_envAdj + treat_SE_oa_2800,
       angle = 90, len = 0.05, code = 3, xpd = NA, lwd = 2,
       col=col_perm[4])

# Significance symbols based on planned comparisons
text(x=-1,y=1.20,label="Trt",cex=1.8,xpd=NA,srt=90)
lines(x=c(0,80),y=c(1,1),lwd=2,xpd=NA)
text(x=-1,y=0.80,label="Env",cex=1.8,xpd=NA,srt=90)
#Planned tukey test with treatment

text(x=10,
     y=1.15,label="***",
     cex = 2,col=col_perm[4], xpd = NA)
text(x=22,
     y=1.15,label="**",
     cex = 2,col=col_perm[4], xpd = NA)
# Planned t tests with environment
print(tTest_df)
text(x=2,
     y=0.90,label="*",
     cex = 2,col=col_perm[2], xpd = NA)
text(x=5,
     y=0.90,label="*",
     cex = 2,col=col_perm[2], xpd = NA)
text(x=10,
     y=0.90,label="*",
     cex = 2,col=col_perm[2], xpd = NA)
text(x=22,
     y=0.90,label="**",
     cex = 2,col=col_perm[2], xpd = NA)
text(x=50,
     y=0.90,label="**",
     cex = 2,col=col_perm[2], xpd = NA)
text(x=80,
     y=0.90,label="**",
     cex = 2,col=col_perm[2], xpd = NA)
text(x=2,
     y=0.80,label="",
     cex = 2,col=col_perm[5], xpd = NA)
text(x=5,
     y=0.80,label="**",
     cex = 2,col=col_perm[5], xpd = NA)
text(x=10,
     y=0.80,label="*",
     cex = 2,col=col_perm[5], xpd = NA)
text(x=22,
     y=0.80,label="**",
     cex = 2,col=col_perm[5], xpd = NA)
text(x=50,
     y=0.80,label="",
     cex = 2,col=col_perm[5], xpd = NA)
text(x=80,
     y=0.80,label="*",
     cex = 2,col=col_perm[5], xpd = NA)
text(x=2,
     y=0.70,label="+",
     cex = 2,col=col_perm[4], xpd = NA)
text(x=5,
     y=0.70,label="+",
     cex = 2,col=col_perm[4], xpd = NA)
text(x=10,
     y=0.70,label="+",
     cex = 2,col=col_perm[4], xpd = NA)
text(x=22,
     y=0.70,label="+",
     cex = 2,col=col_perm[4], xpd = NA)
text(x=50,
     y=0.70,label="+",
     cex = 2,col=col_perm[4], xpd = NA)
text(x=80,
     y=0.70,label="*",
     cex = 2,col=col_perm[4], xpd = NA)

# Letter label
text(x=-8,y=1.3,label="B",cex = 2.2, xpd = NA)
