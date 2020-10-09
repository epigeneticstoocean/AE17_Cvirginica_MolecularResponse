#### Script used to analyze the long term EPF pH and calcification data ####
# Figures from this code : 2

## packages
library(dplyr)
library(plyr)
library(multcompView)
library(ggplot2)
library(cowplot)
library(lmerTest)
library(car)
library(RColorBrewer)
pal <- brewer.pal(n = 12, name = 'Paired')
yellow <- brewer.pal(n=9,name = 'YlOrRd')[4]
col_perm <- c(pal[1:2],pal[5:6],pal[12])
col_perm <- c(pal[1:2],pal[5:6],yellow)
# Located in src Analysis/Phenotype folder, will need to set full working directory or setwd()
setwd("/home/downeyam/Github/AE17_Cvirginica_MolecularResponse")
source("src/Accessory/basicR_functions.R")

#### Data ####
bw <- read.delim("data/Phenotype/AE17_CalcificationComplete.csv",sep=",")
#bw <- read.delim("data/Phenotype/CompletePhenotype.csv",sep=",")
#bw <- bw[bw$CompleteRecord==1,] # Trims data to only samples with bouyant weight data
table(bw$pCO2,bw$Timepoint)
bw$pCO2_fac <- as.factor(bw$pCO2) # Turn treatment into a factor
bw$Timepoint_fac <- as.factor(bw$Timepoint) # Turn time point into a factor
bw$PopOrigin<- as.factor(bw$pop)
bw$TankID<- as.factor(bw$tankID)
bw$EPF_pH <- bw$pHTotal
#bw$EPF_pH <- bw$EPF_pH_Total
bw$EPF_envAdj <- bw$EPF_pH-bw$pH_Total_2W # Calculate delta pH


#### Water Chemistry Data ####
wc <- read.delim("data/water_chem/AE17_WeeklySW_ExperimentalExposureSamples_Simple.csv",sep=",",
                 stringsAsFactors = FALSE)

#### Analysis #####

#### Long term EPF ####
## Full model
EPF_measure <- lmer(EPF_pH ~ pCO2_fac*Timepoint_fac + (1|PopOrigin) + (1|TankID),data=bw)
# Model is too complicated
## No random effects (FINAL MODEL)
EPF_measure_red <- lm(EPF_pH ~ pCO2_fac*Timepoint_fac,data=bw)
summary(aov(EPF_measure_red))
#                       Df Sum Sq Mean Sq F value   Pr(>F)    
#pCO2_fac                2 1.6015  0.8007  22.438 1.29e-06 ***
#Timepoint_fac           1 0.0213  0.0213   0.597    0.446    
#pCO2_fac:Timepoint_fac  2 0.0302  0.0151   0.423    0.659    
#Residuals              29 1.0349  0.0357  
TukeyHSD(aov(EPF_measure_red))$pCO2_fac

# Alternative dropped time and interaction which only stengthens effect of pCO2 
# Simple 
EPF_measure_alt <- lm(EPF_pH ~ pCO2_fac,data=bw)
EPF_measure_final_aov <- aov(EPF_measure_alt)
summary(EPF_measure_final_aov)
            

#### Long term relative EPF ####
## Full model
EPF_measure <- lmer(EPF_envAdj ~ pCO2_fac*Timepoint_fac + (1|PopOrigin) + (1|TankID),data=bw)
# Model is too complicated
## No random effects
EPF_measure_red <- lm(EPF_envAdj  ~ pCO2_fac*Timepoint_fac,data=bw)
summary(aov(EPF_measure_red))
#                       Df Sum Sq Mean Sq F value  Pr(>F)   
#pCO2_fac                2 0.5651 0.28253   7.982 0.00173 **
#Timepoint_fac           1 0.0189 0.01892   0.534 0.47062   
#pCO2_fac:Timepoint_fac  2 0.0412 0.02060   0.582 0.56521   
#Residuals              29 1.0266 0.03540                   
# post hoc comparisons
TukeyHSD(aov(EPF_measure_red))$pCO2_fac

# Alternative dropped time and interaction which only stengthens effect of pCO2 
# Simple 
EPF_measure_alt <- lm(EPF_envAdj  ~ pCO2_fac,data=bw)
EPF_measure_final_aov <- aov(EPF_measure_alt)
summary(EPF_measure_final_aov)


#### Calcification vs Environment ####
## Full model
cal_full <- lmer(PercentChangePerDay_DW2_DW3~ pCO2_Cal + (1|PopOrigin) + (1|TankID),data=bw)
## Reduced model final)
cal_analysis_fixed <- lm(PercentChangePerDay_DW2_DW3~pCO2_Cal,data=bw)
summary(cal_analysis_fixed)

# Full Model
calVspH_full <- lmer(PercentChangePerDay_DW2_DW3 ~ EPF_pH+Timepoint_fac + (1|PopOrigin) + (1|TankID),data=bw)
# Reducded Model (final)
calVspH_simple <- lm(PercentChangePerDay_DW2_DW3 ~ EPF_pH,data=bw)
summary(calVspH_simple)

#### Figure ####

#### Panel A - EPF pH measured end trend ####
# Mean environmental pH
c_mean <- mean(wc$pH_Total[wc$PCO2 == 550])
oa_900_mean <- mean(wc$pH_Total[wc$PCO2 == 1000])
oa_2800_mean <- mean(wc$pH_Total[wc$PCO2 == 2800])

#Statistical Model
tHSD <- TukeyHSD(aov(EPF_pH ~ pCO2_fac,data=bw))

# Create letters of significance
agg <- aggregate(EPF_pH ~ pCO2_fac,data = bw,FUN = max)
agg[,2] <- agg[,2]+0.03
Tukey.levels <- tHSD[['pCO2_fac']][,4]
Tukey.labels <- multcompLetters(Tukey.levels)['Letters']
plot.labels <- names(Tukey.labels[['Letters']])
plot.labels.df <- data.frame(plot.labels,tlevels=Tukey.labels)
labels.df <- merge(plot.labels.df,agg, by.x = 'plot.labels', by.y = 'pCO2_fac', sort = FALSE)
names(labels.df) <- c("pCO2_fac","labels","height")
labels.df$height <- labels.df$height + .05

pA <- ggplot(bw,aes(x=pCO2_fac,y=EPF_pH,fill=pCO2_fac)) + 
  geom_hline(aes(yintercept=c_mean,linetype="Ambient"),colour=col_perm[2],size=.8,show.legend =TRUE) + 
  geom_hline(aes(yintercept=oa_900_mean,linetype="OA 900"),colour=col_perm[5],size=.8) + 
  geom_hline(aes(yintercept=oa_2800_mean,linetype="OA 2800"),colour=col_perm[4],size=.8) +
  geom_boxplot() +
  theme_cowplot() + 
  scale_y_continuous(breaks=c(6.50,7.00,7.50,8.00),limits=c(6.45,8.15),labels=c(" 6.50"," 7.00"," 7.50"," 8.00")) +
  scale_fill_manual(values=c(col_perm[2],col_perm[5],col_perm[4]))+
  labs(x="",y=expression(pH[EPF]~(Total)),fill="") +
  #theme(legend.position = "NULL") +
  geom_text(data = labels.df,
            aes(x = pCO2_fac, y = height, label = labels),
            size=5) +
  scale_linetype_manual(name = "Seawater pH", values = c(2,3,4), 
                      guide = guide_legend(override.aes = list(color = c(col_perm[2],col_perm[5],col_perm[4])))) +
  theme(legend.position="null")
pA


#### Panel B - EPF pH relative measured end trend ####

#Statistical Model
summary(aov(bw$EPF_envAdj ~ bw$pCO2_fac))
tHSD_rel <- TukeyHSD(aov(EPF_envAdj~pCO2_fac,data=bw))

# Creating letters of significance
# Create letters of significance
agg <- aggregate(EPF_envAdj ~ pCO2_fac,data = bw,FUN = max)
agg[,2] <- agg[,2]+0.03
Tukey.levels <- tHSD_rel[['pCO2_fac']][,4]
Tukey.labels <- multcompLetters(Tukey.levels)['Letters']
plot.labels <- names(Tukey.labels[['Letters']])
plot.labels.df <- data.frame(plot.labels,tlevels=Tukey.labels)
labels.df <- merge(plot.labels.df,agg, by.x = 'plot.labels', by.y = 'pCO2_fac', sort = FALSE)
names(labels.df) <- c("pCO2_fac","labels","height")
labels.df$height <- labels.df$height + .05

pB <- ggplot(bw,aes(x=pCO2_fac,y=EPF_envAdj,fill=pCO2_fac)) + 
  #geom_hline(yintercept=0,linetype=3,size=0.8) +
  geom_boxplot() +
  theme_cowplot() + 
  scale_fill_manual(values=c(col_perm[2],col_perm[5],col_perm[4]))+
  labs(x="",y=expression(paste(Delta," pH (Total)")),fill="") +
  theme(legend.position = "NULL") +
  scale_y_continuous(limits=c(-1.1,0.3),breaks=c(-0.9,-0.60,-0.30,0,0.3),labels=c("-0.90","-0.60","-0.30","0.00","0.30")) +
  geom_text(data = labels.df,
            aes(x = pCO2_fac, y = height, label = labels),
            size=5)
pB

#### Panel C - Calcification vs. pCO2 ####
#Mean calcification rate 
mean_cal <- aggregate(PercentChangePerDay_DW2_DW3~pCO2,bw,FUN=mean)
se_cal <- aggregate(PercentChangePerDay_DW2_DW3~pCO2,bw,FUN=se)
se_pco2 <- aggregate(pCO2_Cal~pCO2,bw,FUN=sd)

mean_cal <- data.frame(pCO2=mean_cal$pCO2,Rel_Change=mean_cal$PercentChangePerDay_DW2_DW3,
                       ymin=mean_cal$PercentChangePerDay_DW2_DW3-c(se_cal$PercentChangePerDay_DW2_DW3*1.96),
                       ymax=mean_cal$PercentChangePerDay_DW2_DW3+c(se_cal$PercentChangePerDay_DW2_DW3*1.96),
                       xmin=mean_cal$pCO2-c(se_pco2$pCO2_Cal*1.96),
                       xmax=mean_cal$pCO2+c(se_pco2$pCO2_Cal*1.96))

# Linear model (Calcification explained by pCO2)
out <- summary(lm(PercentChangePerDay_DW2_DW3~pCO2,data=bw))
out
p <- ggplot(mean_cal,aes(x=pCO2,y=Rel_Change,shape=as.factor(pCO2),colour=as.factor(pCO2))) +
  geom_abline(slope = out$coefficients[2,1],intercept = out$coefficients[1,1]) +
  geom_hline(aes(yintercept=0),linetype="dotted") +
  geom_point(aes(size=1.5)) +
  ylim(-0.05,0.05) +
  xlim(280,3200) +
  scale_shape_manual(values=c(16,15,17))+
  scale_colour_manual(values=c(col_perm[2],col_perm[5],col_perm[4]))+
  geom_errorbarh(aes(xmin=xmin, xmax=xmax)) + 
  geom_errorbar(aes(ymin=ymin,ymax=ymax),width=75) 

pC <- p + theme_bw(base_size = 16) + 
  labs(x="pCO2 (uatm)",
       y="Calcification (% change per day)") + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position="none")
pC


#### Panel D - Calcification vs. EPF pH ####
bw$pCO2_name <- "NA"
bw$pCO2_name[bw$pCO2 == unique(bw$pCO2)[1]] <-  "Control"
bw$pCO2_name[bw$pCO2 == unique(bw$pCO2)[2]] <-  "Mod. OA" 
bw$pCO2_name[bw$pCO2 == unique(bw$pCO2)[3]] <-  "High OA"
bw$pCO2_name <- factor(bw$pCO2_name,levels = c("Control", "Mod. OA", "High OA"))

p <- ggplot(bw,aes(x=EPF_pH,y=PercentChangePerDay_DW2_DW3,colour=pCO2_name)) +
  geom_point(size=4,aes(shape=pCO2_name)) + 
  ylim(-0.06,0.05) +
  scale_color_manual(values=c(col_perm[2],col_perm[4],col_perm[5])) +
  scale_shape_manual(values=c(16,17,15)) +
  geom_abline(slope = calVspH_simple$coefficients[2],
              intercept = calVspH_simple$coefficients[1]) 
pD <- p + 
  theme_bw(base_size = 16) + 
  labs(x=expression(pH[EPF]~(Total)),
       y="Calcification (% change per day)",
       colour="",
       shape="") + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.spacing.y = unit(0, "mm"),
        legend.position = c(0.2,.85),
        legend.background = element_rect(linetype = 1, 
                                         size = 0.5, 
                                         colour = 1),
        axis.line = element_line(colour = "black")) 
pD

#### Final plot ####
plot_grid(pA,pB,pC,pD,
          ncol=2,
          labels=c("A","B","C","D"))