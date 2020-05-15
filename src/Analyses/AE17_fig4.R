#### This script is used summarize gene expression response to OA and time ####
 # Figures from this code : 4
 # Supp. Figures from this code : XX

#### Multivariate analysis (PCA,PERMANOVA,DAPC) for gene expression and DNA methylation data ####

## packages
knitr::opts_chunk$set(echo = TRUE)
library(ggplot2)
library(cowplot)
library(edgeR)
library(ashr)
library(vegan)
library(sf)
library(adegenet)
library(ggplot2)
library(WGCNA)
library(dplyr)
library(ape)
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
gc <- gc_qualityweights$E
# Using the log2-cpm count matrix produced via the edgeR and limma normalization steps.
# This means the raw counts were first standardized among individuals to 1mill per ind,
# then transformed with a log2.
# FYI the negative values that are created by log transformation of counts < 1 are 
# not an issue for downstream analysis.
## Given what we saw with individual 17005, I removed it as an outlier from the analysis
gc_sans17005 <- gc[,model$ID != "17005"]
model_sans17005 <- model[model$ID != "17005",]

#### GE PERMANOVA ####
# All individuals
(out_gc <- adonis(t(gc)~Time:Treatment+Time+Treatment+Pop+Lane,data=model,
                  permutations = 5000,method = "manhattan" ))
# With 17005 removed
(out_gc <- adonis(t(gc_sans17005)~Time:Treatment+Time+Treatment+Pop+Lane,data=model_sans17005,
                  permutations = 5000,method = "manhattan" ))
# Slightly different results between data sets, but significance is ultimately the same

#### GE PCA #####
## GENE EXPRESSION PCA 
pca <- prcomp(t(gc))
eigs <- pca$sdev^2
temp_df <- data.frame(x=pca$x[,1],y=pca$x[,2],condition=model$SFV)
P3 <- ggplot(temp_df,aes(x=x,y=y,colour=condition)) + 
  geom_point(size=4) + 
  theme_cowplot() + 
  scale_colour_manual(values = col_perm[c(1,2,3,4)],
                      labels=c('Ambient\n     D9','Ambient\n     D80','OA 2800\n     D9', 'OA 2800\n     D80')) + 
  labs(x=paste0("PC1 (",round(eigs[1] / sum(eigs)*100,1),"%)"),
       y=paste0("PC2 (",round(eigs[2] / sum(eigs)*100,1),"%)"),
       title = "",colour="") +
  theme(plot.title = element_text(hjust = 0.5))
P3

#### Alternative running a PCOA rather than a PCA ####
# gc_dist <- vegdist(t(gc_qualityweights$E),method = "manhattan")
# pcoa_out <- pcoa(gc_dist)
# plot(pcoa_out$values$Cumul_eig~rownames(pcoa_out$values))
# df <- data.frame(pcoa1=pcoa_out$vectors[,1],pcoa2=pcoa_out$vectors[,1],model$SFV)
# ggplot(df,aes(x,y,))

#### Alternative with 17005 (Day 80 ambient individual removed) ####
  pca <- prcomp(t(gc_sans17005))
  pca$x
  eigs <- pca$sdev^2
  temp_df <- data.frame(x=pca$x[,1],y=pca$x[,2],
                        condition=model_sans17005$SFV)
  P3_alt <- ggplot(temp_df,aes(x=x,y=y,colour=condition,shape=condition)) +
    geom_point(size=4) +
    theme_cowplot() +
    scale_colour_manual(values = col_perm[c(1,2,3,4)],
                        labels=c('Ambient\n     D9','Ambient\n     D80',
                                 'OA 2800\n     D9', 'OA 2800\n     D80')) +
    scale_shape_manual(values = c(15,15,17,17),
                        labels=c('Ambient\n     D9','Ambient\n     D80',
                                 'OA 2800\n     D9', 'OA 2800\n     D80')) +
    labs(x=paste0("PC1 (",round(eigs[1] / sum(eigs)*100,1),"%)"),
         y=paste0("PC2 (",round(eigs[2] / sum(eigs)*100,1),"%)"),
         title = "",colour="",shape="") +
    theme(plot.title = element_text(hjust = 0.5))
  P3_alt
  
#### GE DAPC ####
### GE DAPC (2 part) ###
# Creating DF by treatment with first timepoint
gc_comb <- gc
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
whole_meta$coord_transform <- log(whole_meta$coord+abs(min(whole_meta$coord))+1)

P4 <- ggplot(whole_meta,aes(x=coord,fill=SFV)) + 
  geom_density(adjust=2) + xlim(-7,7) +
  scale_fill_manual(values = col_perm[c(1,2,3,4)],
                    labels=c('Ambient\n     D9','Ambient\n     D80',
                             'OA 2800\n     D9', 'OA 2800\n     D80')) +
  theme_cowplot() +
  labs(x="Coordinate",
       y="Density",
       title = "") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_blank())

plot_top <- plot_grid(P3,P4,ncol=2, labels=c("A","B"))

### Alternative without 17005 ###
# Creating DF by treatment with first timepoint
gc_comb <- gc_sans17005
early_time_counts_V1 <- gc_comb[,model_sans17005$Day == 9]
early_time_meta <- model_sans17005[model_sans17005$Day == 9,]
#dapc(t(early_time_counts_V1),early_time_meta$treatment)
dapc_treatment <- dapc(t(early_time_counts_V1),early_time_meta$treatment,n.pca=7,n.da=1)
# PCs = 7, clusters = 1
# Map day 80 samples
early_time_meta$coord<- unlist(dapc_treatment$ind.coord[,1])
late_time_counts <- gc_comb[,model_sans17005$Day == 80]
late_time_meta <- model_sans17005[model_sans17005$Day == 80,]
late_time_meta <- late_time_meta
predict_values <- predict.dapc(dapc_treatment,t(late_time_counts))

late_time_meta$coord <-unlist(predict_values$ind.scores[,1])

whole_meta<- rbind(early_time_meta,late_time_meta)
whole_meta$coord_transform <- log(whole_meta$coord+abs(min(whole_meta$coord))+1)

P4_alt <- ggplot(whole_meta,aes(x=coord,fill=SFV)) + 
  geom_density(adjust=2.5) + xlim(-15,12) + 
  scale_fill_manual(values = col_perm[c(1,2,3,4)],
                    labels=c('Ambient\n     D9','Ambient\n     D80',
                             'OA 2800\n     D9', 'OA 2800\n     D80')) +
  scale_y_continuous(breaks = c(0.01,0.1,0.2,0.8,2),limits = c(0,2),trans="sqrt") + 
  theme_cowplot() + 
  labs(x="Coordinate",
       y="Density (sqrt)",
       title = "") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_blank()) 
P4_alt

#### Final Figure ####
plot_grid(P3_alt,P4_alt,labels=c("A","B"),nrow=2)

# Note this figure was modified afterwards in inkscape to add arrows to indicate the net direction
# of movement along the discriminant function from day 9 to day 80.
