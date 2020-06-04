#### Script used to analyze the long term EPF pH and calcification data ####
# Figures from this code : 2

## packages
library(dplyr)
library(plyr)
library(multcompView)
library(ggplot2)
library(lmerTest)
library(car)
library(RColorBrewer)
pal <- brewer.pal(n = 12, name = 'Paired')
col_perm <- c(pal[1:2],pal[5:6],pal[12])
# Located in src Analysis/Phenotype folder, will need to set full working directory or setwd()
setwd("/home/downeyam/Github/AE17_Cvirginica_MolecularResponse")
source("src/Accessory/basicR_functions.R")

#### Data ####

#### Analysis #####

#### Long term EPF ####
## Full model
EPF_measure <- lmer(EPF_pH ~ pCO2_fac*survival_date + (1|PopOrigin) + (1|shelf/tank),data=cal_epf)
# Model is too complicated
## No random effects
EPF_measure_red <- lm(EPF_pH ~ pCO2_fac*as.factor(survival_date),data=cal_epf)
summary(aov(EPF_measure_red))
# Dropped time and interaction 
# Simple (final model)
EPF_measure_final <- lm(EPF_pH ~ pCO2_fac,data=cal_epf)
EPF_measure_final_aov <- aov(EPF_measure_final)
summary(EPF_measure_final_aov)
# post hoc comparisons
TukeyHSD(EPF_measure_final_aov)
plot(TukeyHSD(EPF_measure_final_aov))

#### Long term relative EPF ####
## Full model
EPF_measure <- lmer(EPF_envAdj ~ pCO2_fac*survival_date + (1|PopOrigin) + (1|shelf/tank),data=cal_epf)
# Model is too complicated
## No random effects
EPF_measure_red <- lm(EPF_envAdj  ~ pCO2_fac*as.factor(survival_date),data=cal_epf)
summary(aov(EPF_measure_red))
# Dropped time and interaction 
# Simple (final model)
EPF_measure_final <- lm(EPF_envAdj  ~ pCO2_fac,data=cal_epf)
EPF_measure_final_aov <- aov(EPF_measure_final)
summary(EPF_measure_final_aov)
# post hoc comparisons
TukeyHSD(EPF_measure_final_aov)
plot(TukeyHSD(EPF_measure_final_aov))
anova(EPF_measure_red,EPF_measure_final)

#### Calcification vs Environment ####
## Full model
cal_full <- lmer(calcification~pCO2_calc + (1|PopOrigin) + (1|shelf/tank),data=cal_red2)
## Reduced model final)
cal_analysis_fixed <- lm(calcification~pCO2_calc,data=cal_red2)
summary(cal_analysis_fixed)

#### Calcification vs EPF pH ####
cal_epf <- inner_join(cal_red2,pheno_red)
# Full Model
calVspH_full <- lmer(calcification ~ EPF_pH+sample_date + (1|PopOrigin) + (1|shelf/tank),data=cal_epf)
# Reducded Model (final)
calVspH_simple <- lm(calcification*100 ~ EPF_pH,data=cal_epf)
summary(calVspH_simple)

#### Figure ####

#### Panel A - EPF pH measured end trend ####
# Mean environmental pH
c_mean <- mean(pheno_red$Tank_pH[pheno_red$pCO2_fac==400])
oa_900_mean <- mean(pheno_red$Tank_pH[pheno_red$pCO2_fac==900])
oa_2800_mean <- mean(pheno_red$Tank_pH[pheno_red$pCO2_fac==2800])

#Statistical Model
tHSD <- TukeyHSD(aov(EPF_pH ~ pCO2_fac,data=cal_epf))

# Create letters of significance
agg <- aggregate(EPF_pH ~ pCO2_fac,data = cal_epf,FUN = max)
agg[,2] <- agg[,2]+0.03
Tukey.levels <- tHSD[['pCO2_fac']][,4]
Tukey.labels <- multcompLetters(Tukey.levels)['Letters']
plot.labels <- names(Tukey.labels[['Letters']])
plot.labels.df <- data.frame(plot.labels,tlevels=Tukey.labels)
labels.df <- merge(plot.labels.df,agg, by.x = 'plot.labels', by.y = 'pCO2_fac', sort = FALSE)
names(labels.df) <- c("pCO2_fac","labels","height")


pA <- ggplot(cal_epf,aes(x=pCO2_fac,y=EPF_pH,fill=pCO2_fac)) + 
  geom_hline(aes(yintercept=c_mean,linetype="Ambient"),colour=col_perm[2],size=.8,show.legend =TRUE) + 
  geom_hline(aes(yintercept=oa_900_mean,linetype="OA 900"),colour=col_perm[5],size=.8) + 
  geom_hline(aes(yintercept=oa_2800_mean,linetype="OA 2800"),colour=col_perm[4],size=.8) +
  geom_boxplot() +
  theme_cowplot() + 
  scale_fill_manual(values=c(col_perm[2],col_perm[5],col_perm[4]))+
  labs(x="",y=expression(pH[EPF]~(NBS)),fill="") +
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
summary(aov(cal_epf$EPF_envAdj ~ cal_epf$pCO2_fac))
tHSD_rel <- TukeyHSD(aov(EPF_envAdj~pCO2_fac,data=cal_epf))

# Creating letters of significance
# Create letters of significance
agg <- aggregate(EPF_envAdj ~ pCO2_fac,data = cal_epf,FUN = max)
agg[,2] <- agg[,2]+0.03
Tukey.levels <- tHSD_rel[['pCO2_fac']][,4]
Tukey.labels <- multcompLetters(Tukey.levels)['Letters']
plot.labels <- names(Tukey.labels[['Letters']])
plot.labels.df <- data.frame(plot.labels,tlevels=Tukey.labels)
labels.df <- merge(plot.labels.df,agg, by.x = 'plot.labels', by.y = 'pCO2_fac', sort = FALSE)
names(labels.df) <- c("pCO2_fac","labels","height")

pB <- ggplot(cal_epf,aes(x=pCO2_fac,y=EPF_envAdj,fill=pCO2_fac)) + 
  geom_hline(yintercept=0,linetype=3,size=0.8) +
  geom_boxplot() +
  theme_cowplot() + 
  scale_fill_manual(values=c(col_perm[2],col_perm[5],col_perm[4]))+
  labs(x="",y=expression(paste(Delta," pH (NBS)")),fill="") +
  theme(legend.position = "NULL") +
  geom_text(data = labels.df,
            aes(x = pCO2_fac, y = height, label = labels),
            size=5)
pB

#### Panel C - Calcification vs. pCO2 ####
#Mean calcification rate 
mean_cal <- aggregate(calcification*100~pCO2,cal_red2,FUN=mean)
se_cal <- aggregate(calcification*100~pCO2,cal_red2,FUN=se)
se_pco2 <- aggregate(pCO2_calc~pCO2,cal_red2,FUN=sd)
mean_cal <- data.frame(pCO2=mean_cal$pCO2,Rel_Change=mean_cal$calcification,
                       ymin=mean_cal$calcification-c(se_cal$calcification*1.96),
                       ymax=mean_cal$calcification+c(se_cal$calcification*1.96),
                       xmin=mean_cal$pCO2-c(se_pco2$pCO2_calc*1.96),
                       xmax=mean_cal$pCO2+c(se_pco2$pCO2_calc*1.96))

# Linear model (Calcification explained by pCO2)
out <- summary(lm(calcification*100~pCO2,data=cal_red2))
out
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
cal_epf$pCO2_name <- "NA"
cal_epf$pCO2_name[cal_epf$pCO2 == unique(cal_epf$pCO2)[1]] <-  "Control"
cal_epf$pCO2_name[cal_epf$pCO2 == unique(cal_epf$pCO2)[2]] <-  "Moderate OA" 
cal_epf$pCO2_name[cal_epf$pCO2 == unique(cal_epf$pCO2)[3]] <-  "High OA"
cal_epf$pCO2_name <- factor(cal_epf$pCO2_name,levels = c("Control", "Moderate OA", "High OA"))

p <- ggplot(cal_epf,aes(x=EPF_pH,y=calcification*100,colour=pCO2_name)) +
  geom_point(size=4,aes(shape=pCO2_name)) + 
  ylim(-0.05,0.05) +
  scale_color_manual(values=c(col_perm[2],col_perm[5],col_perm[4])) +
  scale_shape_manual(values=c(16,15,17)) +
  geom_abline(slope = calVspH_simple$coefficients[2],
              intercept = calVspH_simple$coefficients[1]) 

pD <- p + 
  theme_bw(base_size = 16) + 
  labs(x=expression(pH[EPF]~(NBS)),
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
  