#### Multivariate analysis (PCA,PERMANOVA,DAPC) for gene expression and DNA methylation data ####

## packages
knitr::opts_chunk$set(echo = TRUE)
library(gplots)
library(edgeR)
library(ashr)
library(vegan)
library(sf)
library(adegenet)
library(ggplot2)
#devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

library(RColorBrewer)
pal <- brewer.pal(n = 12, name = 'Paired')
col_perm <- c(pal[1:2],pal[5:6],pal[12])

#### DATA ####
## path
inputDir <- "~/Github/AE17_Cvirginica_MolecularResponse/"

setwd(inputDir)
## Meta Data for the Model
model<-readRDS("data/meta/metadata_20190811.RData")
model_dnam<-model[model$ID != "17099",]
#### Gene Matrix ####
gc_qualityweights<-readRDS("data/Analysis/RNA_gene_postVoomAndNormalization_DGEListObj.RData")
#### DNA Methylation ####
library(data.table)
meth_all_meth <- fread("data/MBDBS_seq/methylKitObj_all_cov5Filtered_united_MethylCCounts.csv")
meth_all_total <- fread("data/MBDBS_seq/methylKitObj_all_cov5Filtered_united_totalCounts.csv")

meth_beta <- meth_all_meth/meth_all_total
beta <- as.matrix(meth_beta)
class(beta) <- "numeric"
#### Analysis ####

#### GE PERMANOVA ####
(out_gc <- adonis(t(gc_qualityweights$E)~Time:Treatment+Time+Treatment+Pop+Lane,data=model,
                  permutations = 5000,method = "manhattan" ))
pairwise.adonis(t(gc_qualityweights$E),factors=model$SFV,p.adjust.m = "fdr",sim.function='vegdist',
                sim.method='manhattan',perm = 5000)

#### GE DAPC (2 part) ####
# Creating DF by treatment with first timepoint
gc_comb <- gc_qualityweights$E
early_time_counts_V1 <- gc_comb[,model$Day == 9]
early_time_meta <- model[model$Day == 9,]
#dapc(t(early_time_counts_V1),early_time_meta$treatment)
dapc_treatment <- dapc(t(early_time_counts_V1),early_time_meta$treatment,n.pca=7,n.da=1)
# PCs = 7, clusters = 1
# Map day 80 samples
early_time_meta$coord<- unlist(dapc_treatment$ind.coord[,1])
late_time_counts <- gc_comb[,model$Day == 80]
late_time_meta <- model[model$Day == 80,]

predict_values <- predict.dapc(dapc_treatment,t(late_time_counts))

late_time_meta$coord <-unlist(predict_values$ind.scores[,1])

whole_meta<- rbind(early_time_meta,late_time_meta)
# Separate out data by time x treatment for final plot
D9_400 <- whole_meta[whole_meta$SFV == "09.400",]
D9_400_Density <- density(D9_400$coord)
D9_2800 <- whole_meta[whole_meta$SFV == "09.2800",]
D9_2800_Density <- density(D9_2800$coord)
D80_400 <- whole_meta[whole_meta$SFV == "80.400",]
D80_400_Density <- density(D80_400$coord)
D80_2800 <- whole_meta[whole_meta$SFV == "80.2800",]
D80_2800_Density <- density(D80_2800$coord)

#### GE DAPC (1 part treatment x time) ####
dapc_SFV_10<-dapc(t(gc_qualityweights$E),model$SFV,n.pca=8,n.da=3)
# PCs = 8, clusters = 3
output <- data.frame(Trt=model$Treatment,Time=model$Time,dapc_SFV_10$ind.coord)
#Basic ggplot similar to supplemental 
(nogrouping_ge <- ggplot(output,
                         aes(x=LD1,
                             y=LD2,fill=as.factor(interaction(Trt,Time)),
                             colour=as.factor(interaction(Trt,Time)),
                             shape=as.factor(interaction(Trt,Time)),
                             size=as.factor(interaction(Trt,Time)))) + 
    geom_point(aes()) + #geom_density(alpha=0.1) + #xlim(-28,28) + 
    labs(title="",
         x="Discriminant function 1",
         y="Discriminant function 2",
         colour="Treatment",
         fill="Treatment") +
    theme_classic() +
    scale_color_manual(values=c(col_perm[1],col_perm[3],col_perm[2],col_perm[4])) +
    scale_fill_manual(values=c(col_perm[1],col_perm[3],col_perm[2],col_perm[4])) +
    scale_shape_manual(values=c(16,16,10,10)) +
    scale_size_manual(values=c(5,5,5,5)) +
    labs(colour="Treatment",size="Treatment",shape="Treatment")) +
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=18))

#### DNA Methylation PERMANOVA ####
(out_dnam <- adonis(t(beta)~Treatment*Time+Pop+Lane,data=model_dnam,
                    permutations = 5000,method="manhattan"))
#### DNA Methylation DAPC (2 part) ####
## Creating initial dapc
beta_comb <- beta
early_time_dnam <- beta_comb[,model_dnam$Day == 9]
early_time_meta_dnam <- model_dnam[model_dnam$Day == 9,]

dapc_treatment <- dapc(t(early_time_dnam),early_time_meta_dnam$treatment,n.pca=8,n.da=1)
# PCs = 8 (accounts for ~ 80% of variation)
# clusters = 1

early_time_meta_dnam$coord<- unlist(dapc_treatment$ind.coord[,1])
## Mapping Day 80
late_time_dnam <- beta_comb[,model_dnam$Day == 80]
late_time_meta_dnam <- model_dnam[model_dnam$Day == 80,]

predict_values <- predict.dapc(dapc_treatment,t(late_time_dnam))

late_time_meta_dnam$coord <-unlist(predict_values$ind.scores[,1])

whole_meta_dnam<- rbind(early_time_meta_dnam,late_time_meta_dnam)
## Subseting data for plot
D9_400_dnam<- whole_meta_dnam[whole_meta_dnam$SFV == "09.400",]
D9_400_Density_dnam <- density(D9_400_dnam$coord)
D9_2800_dnam <- whole_meta_dnam[whole_meta_dnam$SFV == "09.2800",]
D9_2800_Density_dnam <- density(D9_2800_dnam$coord)
D80_400_dnam <- whole_meta_dnam[whole_meta_dnam$SFV == "80.400",]
D80_400_Density_dnam <- density(D80_400_dnam$coord)
D80_2800_dnam <- whole_meta_dnam[whole_meta_dnam$SFV == "80.2800",]
D80_2800_Density_dnam <- density(D80_2800_dnam$coord)

#### DNA Methylation DAPC (1 part treatment x time) ####
dapc_all_dnam<-dapc(t(beta),model_dnam$SFV,n.pca=8,n.da=3)
# PCs = 8, clusters = 3
output <- data.frame(Trt=model_dnam$Treatment,Time=model_dnam$Time,dapc_all_dnam$ind.coord)
## Figure similar to what is found in supplemental
(nogrouping_dnam <- ggplot(output,
                           aes(x=LD1,
                               y=LD2,fill=as.factor(interaction(Trt,Time)),
                               colour=as.factor(interaction(Trt,Time)),
                               shape=as.factor(interaction(Trt,Time)),
                               size=as.factor(interaction(Trt,Time)))) + 
    geom_point(aes()) +  
    labs(title="",
         x="Discriminant function 1",
         y="Discriminant function 2",
         colour="Treatment",
         fill="Treatment") +
    theme_classic() +
    scale_color_manual(values=c(col_perm[1],col_perm[3],col_perm[2],col_perm[4])) +
    scale_fill_manual(values=c(col_perm[1],col_perm[3],col_perm[2],col_perm[4])) +
    scale_shape_manual(values=c(16,16,10,10)) +
    scale_size_manual(values=c(5,5,5,5)) +
    labs(colour="Treatment",size="Treatment",shape="Treatment")) +
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=18)) 

#### Final Figure ####
m <- matrix(c(1,2,3,3,4,5),nrow = 3,ncol = 2,byrow = TRUE)
layout(mat = m,heights = c(0.4,0.04,0.4))
par(mar = c(5,5,3,2))
par(xpd=TRUE)

# Presets for panel letter position
let_ge_pca <- c(-185,105)
let_ge_dapc <- c(-11.8,3.7)
let_dnam_pca <- c(-75,38)
let_dnam_dapc <-c(-5,12.3)

#par(mfrow=c(1,1))
## GENE EXPRESSION PCA 
pca <- prcomp(t(gc_qualityweights$E))
eigs <- pca$sdev^2

color_comb <- c(col_perm[1:4]) # colors for population 
model$colors <- "" 
model$colors[model$SFV == unique(model$SFV)[1]] <-  color_comb[2]
model$colors[model$SFV == unique(model$SFV)[2]] <-  color_comb[1]
model$colors[model$SFV == unique(model$SFV)[3]] <-  color_comb[4]
model$colors[model$SFV == unique(model$SFV)[4]] <-  color_comb[3]
model$pch <- 16
model$pch[model$colors == color_comb[2] | model$colors == color_comb[4]] <- 10

ordiplot(pca,type="n",
         xlab="",
         ylab="",
         cex.axis=1)

mtext(paste0("PC1 (",round(eigs[1] / sum(eigs)*100,1),"% Variance Explained)"), side=1, line=2.8, cex=1)
mtext(paste0("PC2 (",round(eigs[2] / sum(eigs)*100,1),"% Variance Explained)"), side=2, line=2.2, cex=1)

points(pca$x[,2]~pca$x[,1],col=model$colors,pch=model$pch,cex=2.5)

ordispider(pca,model$SFV,col = color_comb,lwd=2.5)

text(x=let_ge_pca[1],y=let_ge_pca[2],labels = "A",cex=2)

## DNAm PCA 
pca <- prcomp(t(beta))
eigs <- pca$sdev^2

color_comb <- c(col_perm[1],col_perm[2],col_perm[3],col_perm[4]) # colors for population 
model_dnam$colors <- "" 
model_dnam$colors[model_dnam$SFV == unique(model_dnam$SFV)[1]] <-  color_comb[2]
model_dnam$colors[model_dnam$SFV == unique(model_dnam$SFV)[2]] <-  color_comb[1]
model_dnam$colors[model_dnam$SFV == unique(model_dnam$SFV)[3]] <-  color_comb[4]
model_dnam$colors[model_dnam$SFV == unique(model_dnam$SFV)[4]] <-  color_comb[3]
model_dnam$pch <- 16
model_dnam$pch[model_dnam$colors == color_comb[2] | model_dnam$colors == color_comb[4]] <- 10

ordiplot(pca,type="n",
         xlab="",
         ylab="",
         cex.axis=1)

points(pca$x[,2]~pca$x[,1],col=model_dnam$colors,pch=model_dnam$pch,cex=2)

ordispider(pca,model_dnam$SFV,col = color_comb,lwd=2.5)

mtext(paste0("PC1 (",round(eigs[1] / sum(eigs)*100,1),"% Variance Explained)"), side=1, line=2.8, cex=1)
mtext(paste0("PC2 (",round(eigs[2] / sum(eigs)*100,1),"% Variance Explained)"), side=2, line=2.2, cex=1)

text(x=let_dnam_pca[1],y=let_dnam_pca[2],labels = "B",cex=2)

## Single legend PCA
par(mar = c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x="center",inset=0,horiz = TRUE,
       legend=c("Day 09 : Ambient","Day 09 : 2800 ","Day 80 : Ambient","Day 80 : 2800"),
       col=c(col_perm[1],col_perm[3],col_perm[2],col_perm[4]),
       pch=c(16,16,10,10),
       cex = 1.5)

## GENE EXPRESSION DAPC
par(mar = c(5,5,3,2))
plot(D9_400_Density, main="",xlim=c(-9,9),ylim=c(0,3.2),
     xlab="",ylab="")
mtext(paste0("Discriminant Function 1"), side=1, line=2.8, cex=1)
mtext(paste0("Density"), side=2, line=2.2, cex=1)
polygon(D9_400_Density, col=alpha(col_perm[1],1), border=col_perm[1])
polygon(D9_2800_Density, col=alpha(col_perm[3],1), border=col_perm[4])
polygon(D80_400_Density, col=alpha(col_perm[2],1), 
        border=col_perm[2],density = 20,cex=100)
polygon(D80_2800_Density, col=alpha(col_perm[4],1),
        border=col_perm[4],density = 20,angle=0,cex=5)
text(x=let_ge_dapc[1],y=let_ge_dapc[2],labels = "C",cex=2)

## DNAm DAPC 
plot(D9_400_Density_dnam, main="",xlim=c(-4,3),ylim=c(0,11),
     xlab="",ylab="")
mtext(paste0("Discriminant Function 1"), side=1, line=2.8, cex=1)
mtext(paste0("Density"), side=2, line=2.2, cex=1)
polygon(D9_400_Density_dnam, col=alpha(col_perm[1],1), border=col_perm[1])
polygon(D9_2800_Density_dnam, col=alpha(col_perm[3],1), border=col_perm[3])
polygon(D80_400_Density_dnam, col=alpha(col_perm[2],1), 
        border=col_perm[2],density = 20,cex=100)
polygon(D80_2800_Density_dnam, col=alpha(col_perm[4],1),
        border=col_perm[4],density = 20,angle=0,cex=5)
text(x=let_dnam_dapc[1],y=let_dnam_dapc[2],labels = "D",cex=2)

## Single legend DAPC
par(mar = c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x="top",inset=0,horiz = TRUE,
       legend=c("Day 09 : Ambient","Day 09 : 2800 ","Day 80 : Ambient","Day 80 : 2800"),
       fill = c(alpha(col_perm[1],1),
                alpha(col_perm[3],1),
                alpha(col_perm[2],1),
                alpha(col_perm[4],1)),
       density = c(1000,1000,50,50),angle = c(0,0,45,0),
       cex = 1.5)
## NOTE: LABELS were modified post production to increase size and legability
