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

#### Data ####
epf_exp <- readRDS("data/Phenotype/AE17_summaryPhenotype_exposure.RData")

#### Analysis ####

# Random effects model using lmer

#### Measured pH vs. treatment and time ####
## Full Model
epfAllTP_full <- lmer(EPF_pH~pCO2_fac*timepoint_fac + (1|PopOrigin) + (1|shelf/tank),data=epf_exp) 
## Reduced (final) model 
epfAllTP_red <- lmer(EPF_pH~pCO2_fac*timepoint_fac + (1|shelf/tank),data=epf_exp)
# LRT to check if reduced model is diff. than full model
anova(epfAllTP_red,epfAllTP_full) # It isnt, so we use the simpler reduced model
# ANOVA
anova(epfAllTP_red)
# Planned Comparisons  
# Create 'k' object with planned comparisons.
# These will include comparisons within each timepoint for
# for either OA treatment (900 or 2800) vs the control.
group <- paste0(pheno_red$pCO2_fac,":",pheno_red$timepoint_fac)
mod_matrix <- model.matrix(epfAllTP_red)
agg_mod_matrix <- aggregate(mod_matrix~group,FUN=mean)
rownames(agg_mod_matrix) <- agg_mod_matrix$group
agg_mod_matrix <- agg_mod_matrix[,-1]
lc2 <- as.matrix(agg_mod_matrix)
k <- rbind("400v900_1"=lc2["900:1",]-lc2["400:1",],
           "400v2800_1"=lc2["2800:1",]-lc2["400:1",],
           "400v900_2"=lc2["900:2",]-lc2["400:2",],
           "400v2800_2"=lc2["2800:2",]-lc2["400:2",],
           "400v900_9"=lc2["900:9",]-lc2["400:9",],
           "400v2800_9"=lc2["2800:9",]-lc2["400:9",],
           "400v900_22"=lc2["900:22",]-lc2["400:22",],
           "400v2800_22"=lc2["2800:22",]-lc2["400:22",],
           "400v900_50"=lc2["900:50",]-lc2["400:50",],
           "400v2800_50"=lc2["2800:50",]-lc2["400:50",],
           "400v900_79"=lc2["900:79",]-lc2["400:79",],
           "400v2800_79"=lc2["2800:79",]-lc2["400:79",]
)
(epfAllTP_posthoc_Correction <- summary(glht(epfAllTP_red,linfct=k),adjusted(type = "fdr")))

#### Relative pH vs. treatment and time ####
## Full Model
epfenvAllTP_full <- lmer(EPF_envAdj~pCO2_fac*timepoint_fac + (1|PopOrigin) + (1|shelf/tank),data=epf_exp) 
## Reduced (final) Model
epfenvAllTP_red <- lmer(EPF_envAdj~pCO2_fac*timepoint_fac + (1|shelf/tank),data=epf_exp)
# LRT to check if reduced model is diff. than full model
anova(epfenvAllTP_red,epfenvAllTP_full) # It isnt, so we use the simpler reduced model
# ANOVA 
anova(epfenvAllTP_red)
# Planned comparisons (using same comparisons when looking at measured pH as response variable)
(epfAllTP_posthoc_Correction <- summary(glht(epfenvAllTP_red,linfct=k),adjusted(type = "fdr")))


#### Figure #######################################################################
# This figure uses base plot() to create a two-panel figure with the EPF timeseries

# Summarize the pH of the environment (featured as lines on figure)
c_mean <- mean(epf_exp$Tank_pH[epf_exp$pCO2_fac==400])
oa_900_mean <- mean(epf_exp$Tank_pH[epf_exp$pCO2_fac==900])
oa_2800_mean <- mean(epf_exp$Tank_pH[epf_exp$pCO2_fac==2800])

# Code for setting the plotting space for the two panels (plus middle section for shared legend)
m <- matrix(c(1,2,3),nrow = 3,ncol = 1,byrow = TRUE)
layout(mat = m,heights = c(0.5,0.05,0.5),widths = c(.1))
par(mar = c(5,5,3,2))

#### Panel A - measured EPF pH ####
### Summarize measured EPF pH for plotting
## Take the means for each timepoint and treatment level
treatSeq_means <- aggregate(EPF_pH~timepoint+pCO2_fac,epf_exp,FUN=mean)
treat_means_ctrl <- treatSeq_means[treatSeq_means$pCO2_fac == "400",]
treat_means_oa_900 <- treatSeq_means[treatSeq_means$pCO2_fac == "900",]
treat_means_oa_2800 <- treatSeq_means[treatSeq_means$pCO2_fac == "2800",]
# Stagger the timepoints for 900 and 2800 treatments so they can be seen better on plot
treat_means_oa_900$timepoint <- treat_means_oa_900$timepoint + 1
treat_means_oa_2800$timepoint <- treat_means_oa_2800$timepoint + 2
# Manually alter timepoint 2 placement along x-axis for visual clarity
treat_means_ctrl$timepoint[2] <- 4
treat_means_oa_900$timepoint[2] <- 5
treat_means_oa_2800$timepoint[2] <- 6
## Take the standard error (SE) for each timepoint and treatment level
treatSeq_SE <- aggregate(EPF_pH~timepoint+pCO2_fac,epf_exp,FUN=se)
treat_SE_ctrl <- treatSeq_SE[treatSeq_SE$pCO2_fac == "400",]$EPF_pH
treat_SE_oa_900 <- treatSeq_SE[treatSeq_SE$pCO2_fac == "900",]$EPF_pH
treat_SE_oa_2800 <- treatSeq_SE[treatSeq_SE$pCO2_fac == "2800",]$EPF_pH

# Base plot with 400
bp <- plot(treat_means_ctrl$EPF_pH~treat_means_ctrl$timepoint,
           ylab=expression(pH[EPF]~(NBS)),xlab="Time (Days)",
           col=col_perm[2],pch=16,cex=2,cex.axis = 1.5,cex.lab=1.5,
           ylim = c(6.5,8.1),xlim=c(0,82),bty="n")

# Mean environment lines underneath other lines    
abline(h=c_mean,col=col_perm[2],lty=2) # control
abline(h=oa_900_mean,col=col_perm[5],lty=2) # 900
abline(h=oa_2800_mean,col=col_perm[4],lty=2) # 2800

# 400 EPF lines
lines(treat_means_ctrl$EPF_pH~treat_means_ctrl$timepoint,
      col="lightblue4")
arrows(x0 = treat_means_ctrl$timepoint,
       x1 = treat_means_ctrl$timepoint,
       y0 = c(treat_means_ctrl$EPF_pH - treat_SE_ctrl),
       y1 = treat_means_ctrl$EPF_pH + treat_SE_ctrl,
       angle = 90, len = 0.05,
       code = 3, xpd = NA, lwd = 2,
       col=col_perm[2])
# 900 EPF lines
points(treat_means_oa_900$EPF_pH~treat_means_oa_900$timepoint,
       col=col_perm[5],pch=15,cex=2)
lines(treat_means_oa_900$EPF_pH~treat_means_oa_900$timepoint,
      col=col_perm[5])
arrows(x0 = treat_means_oa_900$timepoint,
       x1 = treat_means_oa_900$timepoint,
       y0 = c(treat_means_oa_900$EPF_pH - treat_SE_ctrl),
       y1 = treat_means_oa_900$EPF_pH + treat_SE_ctrl,
       angle = 90, len = 0.05,
       code = 3, xpd = NA, lwd = 2,
       col=col_perm[5])
# 2800 EPF lines 
points(treat_means_oa_2800$EPF_pH~treat_means_oa_2800$timepoint,
       col=col_perm[4],pch=17,cex=2)
lines(treat_means_oa_2800$EPF_pH~treat_means_oa_2800$timepoint,
      col=col_perm[4])
arrows(x0 = treat_means_oa_2800$timepoint,
       x1 = treat_means_oa_2800$timepoint,
       y0 = c(treat_means_oa_2800$EPF_pH - treat_SE_oa_2800),
       y1 = treat_means_oa_2800$EPF_pH + treat_SE_oa_2800,
       angle = 90, len = 0.05, code = 3, xpd = NA, lwd = 2,
       col=col_perm[4])

# Based on planned comparisons, significance stars
text(x=5+1,y=8,label="*",cex = 2,col=col_perm[4], xpd = NA)
text(x=51+1,y=8,label="*",cex = 2,col=col_perm[4], xpd = NA)
text(x=80+1,y=8,label="*",cex = 2,col=col_perm[4], xpd = NA)

# Letter labels
text(x=-5,y=8.25,label="A",cex = 2.2, xpd = NA)

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

#### Panel B - diff EPF pH #### 
par(mar = c(5,5,3,2))
## Take the means for each timepoint and treatment level
treatSeq_means <- aggregate(EPF_envAdj~timepoint+pCO2_fac,epf_exp,FUN=mean)
treat_means_ctrl <- treatSeq_means[treatSeq_means$pCO2_fac == "400",]
treat_means_oa_900 <- treatSeq_means[treatSeq_means$pCO2_fac == "900",]
treat_means_oa_2800 <- treatSeq_means[treatSeq_means$pCO2_fac == "2800",]
# Stagger the timepoints for 900 and 2800 treatments so they can be seen better on plot
treat_means_oa_900$timepoint <- treat_means_oa_900$timepoint + 1
treat_means_oa_2800$timepoint <- treat_means_oa_2800$timepoint + 2
# Manually alter timepoint 2 placement along x-axis for visual clarity
treat_means_ctrl$timepoint[2] <- 4
treat_means_oa_900$timepoint[2] <- 5
treat_means_oa_2800$timepoint[2] <- 6
## Take the standard error (se) for each timepoint and treatment level
treatSeq_SE <- aggregate(EPF_envAdj~timepoint+pCO2_fac,epf_exp,FUN=se)
treat_SE_ctrl <- treatSeq_SE[treatSeq_SE$pCO2_fac == "400",]$EPF_envAdj
treat_SE_oa_900 <- treatSeq_SE[treatSeq_SE$pCO2_fac == "900",]$EPF_envAdj
treat_SE_oa_2800 <- treatSeq_SE[treatSeq_SE$pCO2_fac == "2800",]$EPF_envAdj

# Base plot with 400
bp <- plot(treat_means_ctrl$EPF_envAdj~treat_means_ctrl$timepoint,
           ylab=expression(paste(Delta," pH (NBS)")),
           xlab="Time (Days)",pch=16,
           col=col_perm[2],cex=2,cex.axis = 1.5,cex.lab=1.5,
           ylim = c(-1,1.1),xlim=c(0,82),bty="n")
# Mean environment lines underneath other lines    
abline(h=0,col="black",lty=2) # control
# 400 EPF lines
lines(treat_means_ctrl$EPF_envAdj~treat_means_ctrl$timepoint,
      col=col_perm[2])
arrows(x0 = treat_means_ctrl$timepoint,
       x1 = treat_means_ctrl$timepoint,
       y0 = c(treat_means_ctrl$EPF_envAdj - treat_SE_ctrl),
       y1 = treat_means_ctrl$EPF_envAdj + treat_SE_ctrl,
       angle = 90, len = 0.05,
       code = 3, xpd = NA, lwd = 2,
       col=col_perm[2])
# 900 EPF lines
points(treat_means_oa_900$EPF_envAdj~treat_means_oa_900$timepoint,
       col=col_perm[5],pch=15,cex=2)
lines(treat_means_oa_900$EPF_envAdj~treat_means_oa_900$timepoint,
      col=col_perm[5])
arrows(x0 = treat_means_oa_900$timepoint,
       x1 = treat_means_oa_900$timepoint,
       y0 = c(treat_means_oa_900$EPF_envAdj - treat_SE_oa_900),
       y1 = treat_means_oa_900$EPF_envAdj + treat_SE_oa_900,
       angle = 90, len = 0.05,
       code = 3, xpd = NA, lwd = 2,
       col=col_perm[5])
# 2800 EPF lines 
points(treat_means_oa_2800$EPF_envAdj~treat_means_oa_2800$timepoint,
       col=col_perm[4],pch=17,cex=2)
lines(treat_means_oa_2800$EPF_envAdj~treat_means_oa_2800$timepoint,
      col=col_perm[4])
arrows(x0 = treat_means_oa_2800$timepoint,
       x1 = treat_means_oa_2800$timepoint,
       y0 = c(treat_means_oa_2800$EPF_envAdj - treat_SE_oa_2800),
       y1 = treat_means_oa_2800$EPF_envAdj + treat_SE_oa_2800,
       angle = 90, len = 0.05, code = 3, xpd = NA, lwd = 2,
       col=col_perm[4])

# Significance symbols based on planned comparisons
text(x=-1,y=1.10,label="Trt",cex=1.8,xpd=NA,srt=90)
lines(x=c(0,80),y=c(0.95,0.95),lwd=2,xpd=NA)
text(x=-1,y=0.80,label="Env",cex=1.8,xpd=NA,srt=90)
#Planned tukey test with treatment
text(x=2,
     y=1.05,label="**",
     cex = 2,col=col_perm[4], xpd = NA)
text(x=5,
     y=1.05,label="+",
     cex = 2,col=col_perm[4], xpd = NA)
text(x=10,
     y=1.05,label="***",
     cex = 2,col=col_perm[4], xpd = NA)
text(x=23,
     y=1.05,label="**",
     cex = 2,col=col_perm[4], xpd = NA)
text(x=50,
     y=1.05,label="+",
     cex = 2,col=col_perm[4], xpd = NA)
text(x=2,
     y=1.15,label="*",
     cex = 2,col=col_perm[5], xpd = NA)
text(x=10,
     y=1.15,label="+",
     cex = 2,col=col_perm[5], xpd = NA)
text(x=50,
     y=1.15,label="+",
     cex = 2,col=col_perm[5], xpd = NA)
# Planned t tests with environment
text(x=2,
     y=0.90,label="*",
     cex = 2,col=col_perm[2], xpd = NA)
text(x=5,
     y=0.90,label="+",
     cex = 2,col=col_perm[2], xpd = NA)
text(x=10,
     y=0.90,label="***",
     cex = 2,col=col_perm[2], xpd = NA)
text(x=23,
     y=0.90,label="***",
     cex = 2,col=col_perm[2], xpd = NA)
text(x=50,
     y=0.90,label="*",
     cex = 2,col=col_perm[2], xpd = NA)
text(x=23,
     y=0.80,label="+",
     cex = 2,col=col_perm[5], xpd = NA)
text(x=10,
     y=0.70,label="***",
     cex = 2,col=col_perm[4], xpd = NA)

# Letter label
text(x=-5,y=1.3,label="B",cex = 2.2, xpd = NA)
