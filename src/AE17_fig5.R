#### Gene Expression and DNA Methylation Day 9 ambient comparison####
## Used to compare general relationship between DNA methylation and gene expression

## package
knitr::opts_chunk$set(echo = TRUE)
library(rtracklayer) # on bioconducter
library(dplyr)
library(edgeR)
library(factoextra)
library(cowplot)
library(wesanderson)
library(viridis)
library(adegenet)
library(ggplot2)
library(matrixStats)
library(grid)
library(gridExtra)

library(RColorBrewer)
pal <- brewer.pal(n = 3, name = "PRGn")

#### Functions ####

# Function for subset 'mat' for the target expression,DNA methylation, and gene attributes
subMat_control <- function(x){
  y <- data.frame(gene_length=log2(x$gene_length),
                  exon=log2(x$exon_count),totalCpG=log2(x$all_count),
                  Methylation=x$Mean_9C,Methylation_CV=x$cv_9C,
                  Gene_Expression=x$gene_mean_9C,Gene_Expression_CV=x$gene_cv_9C)
  return(y)
}

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

perCentCpGCoverage <- function(x){sum(c(mat$cov5_count/mat$all_count)*100 > x)}

# Function for generating three panel contribution plot for the different PCs
suppContriFig <- function(x){
  temp <- list()
  for(i in 1:3){
    temp[[i]] <- fviz_contrib(x, choice = "var", axes = i, top = 10,color = "black",fill = "grey",
                              sort.val = "none")
  }
  return(temp)
}

#### Data ####
## path
inputDir <- "~/Github/AE17_Cvirginica_MolecularResponse"
outputDir <- ""

setwd(inputDir)

# Read in meta data 
meta_original <- readRDS("data/meta/metadata_20190811.RData")
meta <- meta_original[meta_original$ID != "17099",]

# Read in Gene Expression Matrix
ge<-readRDS("data/Analysis/RNA_gene_postVoomAndNormalization_DGEListObj.RData")
ge_counts <- ge$E # Extract just the counts
ge_counts <- ge_counts[,colnames(ge_counts) != "17099"] # Remove problematic individual 

# Read in DNA Methylation Data for Gene Bodys
dnam <- readRDS("data/Analysis/DNAm_20200202_AllCountsList_cov5_byFeature.RData")
dnam <-dnam$beta$gene # Only going to look at the beta (dna methylation value) for genes
dnam <- as.matrix(dnam[,2:ncol(dnam)]) # Change data into matrix 
dnam[is.na(dnam)] <- 0
class(dnam) <- "numeric"

# Gene Atribute, Expression and Methylation Summary Table
mat <- readRDS("data/Analysis/Multi_geneSummaryReduced.RData")
# We have overlapping expression and DNA Methylation data for 9626 genes

#### Subset data by CpG Coverage ####

# 20% CpG coverage within a gene
mat_20 <-  mat[mat$cov5_count/mat$all_count>=0.2,]
# 50% CpG coverage within a gene
mat_50 <- mat[mat$cov5_count/mat$all_count>=0.5,]
# 80% CpG coverage within a gene
mat_80 <- mat[mat$cov5_count/mat$all_count>=0.8,]

#The counts should be log transformed
histogram(log2(mat$exon_count))
histogram(log2(mat$gene_length))
histogram(log2(mat$all_count))

## Creates matrix with summary stats for samples collected on day 9 at the control 
# treatment. This is to look at basic DNA methylation - Gene expression correlations
mat_0_control <- subMat_control(mat)
mat_20_control <- subMat_control(mat_20)
mat_50_control <- subMat_control(mat_50)
mat_80_control <- subMat_control(mat_80)

#### Gene Attribut PCA ####
### Control Timepoint (developmental response)
# We are using log2 transformed and scaled variables (except methylation is not log2 transformed)
pca_C <- scale(mat_20_control)
# Perform pca
pca_C_obj <- prcomp(pca_C)
# Scree plot 
fviz_eig(pca_C_obj,addlabels = TRUE)

### Among treatment (ecological response)
pca_amongTrt <- scale(mat_20_amongTrt)
pca_obj_amongTrt <- prcomp(pca_amongTrt)
fviz_eig(pca_obj_amongTrt,addlabels = TRUE)
## Timepoint 9
## Timepoint 80
# pca <- scale(mat_20_control)
# pca_obj <- prcomp(pca)
# fviz_eig(pca_obj,addlabels = TRUE)

#### PCA Attribute Figures ####
## Control
p1_C<-fviz_pca_var(pca_C_obj,col.circle = "white",
                          title="",
                          repel = TRUE # Avoid text overlapping
) + theme_cowplot() +
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=15))
pPCA_C <- plot_grid(p1_C,labels="A")  
pPCA_C

dats <- c("mat","mat_20","mat_50","mat_80")
namePlot <- c("A","B","C ","D")

# Mean Methylation in our control groups (Day 9) vs Gene Expression
pL_mean <- list()
for(i in 1:length(dats)){
  df <- get(dats[i])
  df$density <- get_density(df$gene_mean_9C, df$Mean_9C*100, n = 100)
  temp <- ggplot(df,aes(x=gene_mean_9C,y=Mean_9C*100)) + 
    geom_point(aes(color = density)) + ylim(0,100) + 
    geom_smooth(method=lm,colour="orange",size=2,na.rm = TRUE) +
    scale_color_viridis() +
    labs(x="Gene Expression (log2-cpm)",y="DNA Methylation (%)") +
    theme_cowplot()
  pL_mean[[i]]<- temp
}
# Mean Methylation in our control groups (Day 9) vs Gene Expression CV (log)
pL_cv <- list()
for(i in 1:length(dats)){
  df <- get(dats[i])
  df$density <- get_density(log(df$gene_cv_9C), df$Mean_9C*100, n = 100)
  temp <- ggplot(df,aes(x=log(gene_cv_9C),y=Mean_9C*100)) + 
    geom_point(aes(color = density)) + ylim(0,100) + scale_color_viridis() +
    geom_smooth(method=lm,colour="orange",size=2) +
    labs(x=bquote("Gene Expression"~CV[Ind]~"(log)"),y="DNA Methylation (%)") +
    theme_cowplot()
  pL_cv[[i]]<- temp
}


# Development 
right <- plot_grid(pL_mean[[2]],pL_cv[[2]],labels=c("B","C"),ncol=1)
plot_grid(pPCA_C,right)
ggsave(plot_grid(pPCA_C,right),filename="results/figures/Figure5/Figure5.pdf")
ggsave(plot_grid(pPCA_C,right),filename="results/figures/Figure5/Figure5.png")

## Model
# Mean Methylation vs. Mean Gene Expression
summary(lm(mat_20$gene_mean~mat_20$Mean_9C))
# Mean Methylation vs. CV Gene Expression
summary(lm(mat_20$gene_cv_9C~mat_20$Mean_9C))

length(mat_20$Mean_9C)

#### Supplemental Plots for figure 5 ####

## Percent CpG cover within genes X number of genes ##
# Lets check out different filtering levels based on % CpGs covered in a gene
covRange <- data.frame(Percent_Coverage=seq(0,100,5),
                       NumGenes=sapply(seq(0,100,5),perCentCpGCoverage))
# Plot shows number of genes in dataset vs. % CpGs covered within the gene
coverage_fig <- ggplot(covRange,aes(y=NumGenes,x=Percent_Coverage,colour=viridis(1))) +
  geom_vline(xintercept = 20,colour="purple",size=2) + 
  geom_line(colour="black") + geom_point(size=5) + 
  theme_cowplot() + 
  labs(x="Percent Coverage",y="Number of Genes",
       title="Within Gene CpG Coverage") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))
coverage_fig
ggsave(plot_grid(coverage_fig),filename="results/figures/Figure5/supp_CpGByGeneCoverage.pdf")
ggsave(coverage_fig,filename="results/figures/Figure5/supp_CpGByGeneCoverage.png")

## PC axes % contribution figure
# Control
contribute_C_fig <- suppContriFig(pca_C_obj) 
plot_grid(plotlist =contribute_C_fig, labels = c('A', 'B','C'), label_size = 12,nrow = 3)
ggsave(plot_grid(coverage_fig),filename="results/figures/Figure5/supp_PCAAttContr.pdf")
ggsave(coverage_fig,filename="results/figures/Figure5/supp_PCAAttContr.png")

## Gene Expression
plot_grid(plotlist=pL_mean,labels = namePlot)
## CV
plot_grid(plotlist=pL_cv,labels = namePlot)
