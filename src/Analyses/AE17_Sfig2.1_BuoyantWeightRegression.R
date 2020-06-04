#### Bouyant vs Dry Weight Regression ####

#### Data ####
setwd("/home/downeyam/Github/AE17_Cvirginica_MolecularResponse/data/Phenotype/")
bw <- read.csv("AE17_Buoyant_Dry_Ratio.csv")

#### Analysis and plot ####
plot(Dry~Buoyant,data=bw,xlim=c(0,70),ylim=c(0,130),xlab=c("Buoyant Weight (mg)"),ylab=c("Dry Weight (mg)"))
lm.out <- lm(Dry~Buoyant,data=bw)
(lm.sum <- summary(lm.out))
abline(lm.out)
text(x = 35,y = 35,paste0("y = ",round(lm.sum$coef[2,1],2),"x",round(lm.sum$coef[1,1],2)))
text(x = 35,y = 27,bquote(R^2 == .(round(lm.sum$adj.r.squared,2))))
