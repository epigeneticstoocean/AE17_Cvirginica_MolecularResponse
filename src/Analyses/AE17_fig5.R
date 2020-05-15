#### Gene Expression and DNA Methylation Day 9 ambient comparison ####
## Used to compare general relationship between DNA methylation and gene expression
# Figures from this code : 4
# Supp. Figures from this code : XX

## package
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
hist(log2(mat$exon_count))
hist(log2(mat$gene_length))
hist(log2(mat$all_count))

## Creates matrix with summary stats for samples collected on day 9 at the control 
# treatment. This is to look at basic DNA methylation - Gene expression correlations
mat_0_control <- subMat_control(mat)
mat_20_control <- subMat_control(mat_20)
mat_50_control <- subMat_control(mat_50)
mat_80_control <- subMat_control(mat_80)
mats_list <- list(mat_0_control,mat_20_control,mat_50_control,mat_80_control)

# Labels (used later)
dats <- c("mat","mat_20","mat_50","mat_80")
dats_V2 <- c("No filter","20% coverage","50% coverage","80% coverage")
namePlot <- c("A","B","C ","D")

#### Gene Attribut PCA ####
### Control Timepoint (developmental response)
# We are using log2 transformed and scaled variables (except methylation is not log2 transformed)

# prcomp list
prcomp_obj <- list()
# scree plot list
screeplot_obj <- list()
# pca object
pca_obj <- list()

for(i in 1:length(mats_list)){
# Perform pca
  prcomp_obj[[i]] <- prcomp(scale(mats_list[[i]]))
# Scree plot 
  screeplot_obj[[i]] <- fviz_eig(prcomp_obj[[i]],addlabels = TRUE)
## PCA plot
  pca_obj[[i]] <- fviz_pca_var(prcomp_obj[[i]],col.circle = "white",
                          title="",axes=c(1,2),
                          repel = TRUE # Avoid text overlapping
                          ) + 
    theme_cowplot() +
    theme(axis.title = element_text(size=15),
          axis.text = element_text(size=15))
}
pca_obj[[1]]

## Decided on the 20% coverage matrix [[2]]
(pPCA_C <- plot_grid(pca_obj[[2]],labels="A"))

pca_obj_PC2and3 <- fviz_pca_var(prcomp_obj[[2]],col.circle = "white",
                             title="",axes=c(1,3),
                             repel = TRUE # Avoid text overlapping
) + 
  theme_cowplot() +
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=15))
pca_obj_PC2and3

pPCA_alt <- plot_grid(pca_obj[[2]],pca_obj_PC2and3,labels=c("A","B"),nrow=2)

#### Linear correlations between gene expression and DNA methylation ####
prin_PCA <- princomp(scale(mats_list[[2]]))
pr_PCA <- prcomp(scale(mats_list[[2]]))
prin_PCA$loadings
pr_PCA$rotation

get_eigenvalue(prcomp_obj[[2]])
get_pca(prcomp_obj[[2]])$contrib
?get_pca()
# % contribution calculat as the normalized absolute coordinate value * the square root of cos2
cd2 <- get_pca(prcomp_obj[[2]])$coord[,1]
cs2 <- get_pca(prcomp_obj[[2]])$cos2[,1]
(abs(cd2)*sqrt(cs2))/sum((abs(cd2)*sqrt(cs2)))*100
?prcomp()
prcomp_out <- prcomp_obj[[2]]
prcomp_out$rotation
(prcomp_out$sdev)^2
summary(prcomp_out)
loadings(prcomp_out)
scaleV <- 1#sqrt(prcomp_obj[[2]]$scale)
eigen_value <- abs(prcomp_obj[[2]]$rotation[,1])
eigen_vector <- 1 #prcomp_obj[[2]]$sdev
((eigen_value*scaleV)/sum(eigen_value*scaleV))*100
?fviz_contrib
fviz_contrib(prcomp_obj[[2]], choice = "var", axes = 1, top = 10,color = "black",fill = "grey",
             sort.val = "none")
temp <- suppContriFig(prcomp_obj[[2]])
plot_grid(plotlist =temp,labels=c(namePlot[i]),nrow = 1)

## Mean Methylation in our control groups (Day 9) vs Gene Expression
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
## Mean Methylation in our control groups (Day 9) vs Gene Expression CV (log)
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
## Figure 5 alt
right_alt <- plot_grid(pL_mean[[2]],pL_cv[[2]],labels=c("C","D"),ncol=1)
ggsave(plot_grid(pPCA_alt,right_alt),filename="results/figures/Figure5/Figure5_alt.pdf")
ggsave(plot_grid(pPCA_alt,right_alt),filename="results/figures/Figure5/Figure5_alt.png")
### Statistics ###
  # Linear model with Mean gene DNA methylation as predictor 
  # and Mean gene expression or CV as response  

# Mean Methylation vs. Mean Gene Expression
(cpg_ge <- summary(lm(mat_20$gene_mean~mat_20$Mean_9C)))
# Mean Methylation vs. CV Gene Expression
(cpg_cv <- summary(lm(mat_20$gene_cv_9C~mat_20$Mean_9C)))

# Number of genes included
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
ggsave(plot_grid(coverage_fig),filename="results/figures/Figure5/supp_Fig5_CpGByGeneCoverage.pdf")
ggsave(coverage_fig,filename="results/figures/Figure5/supp_Fig5_CpGByGeneCoverage.png")
  
## Additional PCA figures
suppPCA_fig_list <- list()

for(i in 1:length(pca_obj)){
  suppPCA_fig_list[[i]] <- plot_grid(pca_obj[[i]],labels=c(namePlot[i]),nrow = 1)
}

suppPCA_fig_rows <- list()
for(i in 1:length(dats_V2)){
  temp <- textGrob(dats_V2[i],
                   gp=gpar(col="black",
                           fontsize=15))
  suppPCA_fig_rows[[i]] <-grid.arrange(arrangeGrob(suppPCA_fig_list[[i]],top = temp))
}
suppPCA_final <- plot_grid(plotlist = suppPCA_fig_rows,nrow=2)
ggsave(suppPCA_final,filename="results/figures/Figure5/supp_Fig5_PCAsupp.pdf")
ggsave(suppPCA_final,filename="results/figures/Figure5/supp_Fig5_PCAsupp.png")

## PC axes % contribution figure ##
#Contribution plot list for different coverage thresholds
contribute_fig_list <- list()

for(i in 1:length(prcomp_obj)){
  temp <- suppContriFig(prcomp_obj[[i]]) 
  contribute_fig_list[[i]] <- plot_grid(plotlist =temp,labels=c(namePlot[i]),nrow = 1)
}

contribute_fig_rows <- list()
for(i in 1:length(dats_V2)){
  temp <- textGrob(dats_V2[i],
                   gp=gpar(col="black",
                           fontsize=15))
  contribute_fig_rows[[i]] <-grid.arrange(arrangeGrob(contribute_fig_list[[i]],top = temp))
}
contribute_final <- plot_grid(plotlist = contribute_fig_rows,nrow=length(dats_V2))
ggsave(contribute_final,filename="results/figures/Figure5/supp_Fig5_PCAAttContr.pdf")
ggsave(contribute_final,filename="results/figures/Figure5/supp_Fig5_PCAAttContr.png")

## Gene Expression
plot_grid(plotlist=pL_mean,labels = namePlot)
## CV
plot_grid(plotlist=pL_cv,labels = namePlot)
