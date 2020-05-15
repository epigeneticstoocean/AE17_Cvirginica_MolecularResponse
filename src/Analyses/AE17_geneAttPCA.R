#### Gene Expression, Methylation, and Gene Attribute PCA Analysis ####

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

# Read in list of significant DMLs
dmls <- read.csv("data/Analysis/DNAm_methylkit_dm50_q01_sigSummaryTable.csv")
dmls <- dmls[dmls$X1 == "DML",]
dmls$label <- paste0(dmls$chr,"_",dmls$gene_str,"_",dmls$gene_end) 
dmls_9 <- dmls[dmls$X2 == "tp9",c(35,10,17,30)]
dmls_9 <- dmls_9[!is.na(dmls_9$loc),]
dmls_80 <- dmls[dmls$X2 == "tp80",c(35,10,17,30)]
dmls_80 <- dmls_80[!is.na(dmls_80$loc),]
mat_red <- subset(mat,select=c("label","diff_tp9_Trt","diff_tp80_Trt",
                               "gene_diff_tp9_Trt","gene_diff_tp80_Trt",
                               "gene_cv_9C","gene_cv_9E","gene_cv_80C","gene_cv_80E"))

dmls_9_join <- left_join(dmls_9,mat_red,by="label")
dmls_80_join <- left_join(dmls_80,mat_red,by="label")

plot(dmls_9_join$meth.diff~dmls_9_join$gene_diff_tp9_Trt)
plot(dmls_80_join$meth.diff~dmls_80_join$gene_diff_tp80_Trt)

plot(dmls_9_join$meth.diff~c(dmls_9_join$gene_cv_9E-dmls_9_join$gene_cv_9C))
plot(dmls_80_join$meth.diff~c(dmls_80_join$gene_cv_80E-dmls_80_join$gene_cv_80C))

dmls_9_join$diff.cv <- c(dmls_9_join$gene_cv_9E-dmls_9_join$gene_cv_9C)
dmls_80_join$diff.cv <- c(dmls_80_join$gene_cv_80E-dmls_80_join$gene_cv_80C)

dmls_9_join$methylDirection <- ifelse(dmls_9_join$meth.diff < 0,
                                      "Hypomethylation",
                                      "Hypermethylation")
dmls_80_join$methylDirection <- ifelse(dmls_80_join$meth.diff < 0,
                                       "Hypomethylation",
                                       "Hypermethylation")

p1<-ggplot(dmls_9_join,
           aes(y=meth.diff,x=gene_diff_tp9_Trt,colour=methylDirection)) + 
  geom_point(size=3) + ylim(-100,100) + 
  theme_cowplot() + scale_color_manual(values=pal[c(3,1)]) +
  labs(x=bquote("Gene expression (FoldChange"~Log[2]~")"),y="Methylation difference of DMLs",title="Day 9") +
  theme(plot.title = element_text(hjust=0.5))
p1_nolab <- p1 + theme(axis.title = element_blank(),
                       legend.position = "none")
p1_nolab
p2<-ggplot(dmls_80_join,
           aes(y=meth.diff,x=gene_diff_tp80_Trt,colour=methylDirection)) + 
  geom_point(size=3) + ylim(-100,100) +
  theme_cowplot() + scale_color_manual(values=pal[c(3,1)]) +
  labs(x=bquote("Gene expression (FoldChange"~Log[2]~")"),y="Methylation difference of DMLs",title="Day 80") +
  theme(plot.title = element_text(hjust=0.5))
p2_nolab <- p2 + theme(axis.title = element_blank(),
                       legend.position = "none")
p2_nolab
p3<-ggplot(dmls_9_join,
           aes(y=meth.diff,x=diff.cv,colour=methylDirection)) + 
  geom_point(size=3) + ylim(-100,100) + 
  theme_cowplot() + scale_color_manual(values=pal[c(3,1)]) +
  labs(x="Log2 gene expression CV difference",y="Methylation difference of DMLs",title="Day 9") +
  theme(plot.title = element_text(hjust=0.5))
p3_nolab <- p3 + theme(axis.title = element_blank(),
                       legend.position = "none")
p3_nolab
p4<-ggplot(dmls_80_join,
           aes(y=meth.diff,x=diff.cv,colour=methylDirection)) + 
  geom_point(size=3) + ylim(-100,100) + 
  theme_cowplot() + scale_color_manual(values=pal[c(3,1)]) +
  labs(x="Log2 gene expression CV difference",y="Methylation difference of DMLs",title="Day 80") +
  theme(plot.title = element_text(hjust=0.5))
p4_nolab <- p4 + theme(axis.title = element_blank(),
                       legend.position = "none")
p4_nolab

plot1 <- plot_grid(p1_nolab,p2_nolab,ncol=2)
plot2 <- plot_grid(p3_nolab,p4_nolab,ncol=2)
y.grob1 <- textGrob("Methylation difference of DMLs (%)", 
                   gp=gpar(col="black", fontsize=15), rot=90)

x.grob1 <- textGrob("Gene expression difference (log2 fold change)", 
                   gp=gpar( col="black", fontsize=15))
x.grob2 <- textGrob("Gene expression CV difference (log2)", 
                    gp=gpar( col="black", fontsize=15))
save1<-grid.arrange(arrangeGrob(plot1,bottom = x.grob1))
save2<-grid.arrange(arrangeGrob(plot2, bottom = x.grob2))
plot3<-plot_grid(save1,save2,nrow=2)
saveG<-grid.arrange(arrangeGrob(plot3, left = y.grob1))
saveG

## Cadherin
# target <- "LOC111106820"
# exp_target <- ge$E[rownames(ge$E) == target,]
# boxplot(exp_target~meta_original$SFV)

## Basic Diff Methylation plots vs expression plots among treatments
diff_methylation <- ((mat_50$Mean_9E-mat_50$Mean_9C))
gene_diff <- (mat_50$gene_mean_9E-mat_50$gene_mean_9C)
plot(diff_methylation~gene_diff)
summary(lm(diff_methylation~gene_diff))

diff_methylation <- ((mat_50$Mean_80E-mat_50$Mean_80C))
gene_diff <- (mat_50$gene_mean_80E-mat_50$gene_mean_80C)
plot(diff_methylation~gene_diff)
summary(lm(diff_methylation~gene_diff))

diff_methylation <- (mat_20$Mean)
gene_diff <- log(mat_20$gene_cv_mean_AmongTrt)
plot(diff_methylation~gene_diff)
summary(lm(diff_methylation~gene_diff))

#### Additional filtering #### 
# Lets check out different filtering levels based on % CpGs covered in a gene
perCentCpGCoverage <- function(x){sum(c(mat$cov5_count/mat$all_count)*100 > x)}
covRange <- data.frame(Percent_Coverage=seq(0,100,5),
                       NumGenes=sapply(seq(0,100,5),perCentCpGCoverage))


# Plot shows number of genes in dataset vs. % CpGs covered within the gene
coverage_fig <- ggplot(covRange,aes(y=NumGenes,x=Percent_Coverage,colour=viridis(1))) +
  geom_line(colour="black") + geom_point(size=5) + 
  theme_cowplot() + 
  labs(x="Percent Coverage",y="Number of Genes",
       title="Within Gene CpG Coverage") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))
coverage_fig
ggsave(plot_grid(coverage),filename="results/figures/supp_CpGByGeneCoverage.pdf")
ggsave(coverage,filename="results/figures/supp_CpGByGeneCoverage.png")



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

# Function for subset 'mat' for the target expression,DNA methylation, and gene attributes
subMat_control <- function(x){
  y <- data.frame(gene_length=log2(x$gene_length),
             exon=log2(x$exon_count),totalCpG=log2(x$all_count),
             Methylation=x$Mean_9C,Methylation_CV=x$cv_9C,
             Gene_Expression=x$gene_mean_9C,Gene_Expression_CV=x$gene_cv_9C)
  return(y)
}

subMat_amongTrt_9 <- function(x){
  y <- data.frame(gene_length=log2(x$gene_length),
                  exon=log2(x$exon_count),totalCpG=log2(x$all_count),
                  Diff.Meth=x$diff_tp9_Trt,Methylation_CV=log2(x$cv_mean_AmongTrt),
                  Diff.Exp=x$gene_diff_tp9_Trt,Gene_Expression_CV=log2(x$gene_cv_mean_AmongTrt))
  return(y)
}
subMat_amongTrt_80 <- function(x){
  y <- data.frame(#gene_length=log2(x$gene_length),
                  #exon=log2(x$exon_count),totalCpG=log2(x$all_count),
                  Diff.Meth=x$diff_tp9_Trt,Methylation_CV=log2(x$cv_mean_AmongTrt),
                  Gene_Expression=x$gene_mean,Gene_Expression_CV=log2(x$gene_cv_mean_AmongTrt))
  return(y)
}
x <- mat

plot(mat$gene_cv_mean_AmongTrt~mat$gene_cv_9C)
hist(mat$gene_mean)
## Creates matrix with summary stats for samples collected on day 9 at the control 
 # treatment. This is to look at basic DNA methylation - Gene expression correlations
mat_0_control <- subMat_control(mat)
mat_20_control <- subMat_control(mat_20)
mat_50_control <- subMat_control(mat_50)
mat_80_control <- subMat_control(mat_80)
## Creats matrix with the summary stats for all samples and variaation among treatments
 # to see if there are correlated shifts in DNA methylatoin and gene expression.
mat_0_amongTrt <- subMat_amongTrt_9(mat)
mat_20_amongTrt <- subMat_amongTrt_9(mat_20)
head(mat_20_amongTrt)
mat_50_amongTrt <- subMat_amongTrt_9(mat_50)
mat_80_amongTrt <- subMat_amongTrt_9(mat_80)

#### PCA ####

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
p1_amongTrt<-fviz_pca_var(pca_obj_amongTrt,col.circle = "white",
                 title="",
                 repel = TRUE # Avoid text overlapping
) + theme_cowplot() +
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=15))
# p2_amongTrt <- fviz_pca_ind(p1_amongTrt,title="",xlab="",ylab="",
#                    label = "none") + # hide individual labels 
#   theme_cowplot() + 
#   theme(axis.text.y   = element_blank(),
#         axis.text.x   = element_blank(),
#         panel.border = element_rect(colour = "black", fill=NA, size=1))

pPCA_amongTrt <- plot_grid(p1_amongTrt,labels="A")  
pPCA_amongTrt

## Among Treatment
p1<-fviz_pca_var(pca_C_obj,col.circle = "white",
                 title="",
                 repel = TRUE # Avoid text overlapping
) + theme_cowplot() +
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=15))
# p2 <- fviz_pca_ind(pca_C_obj,title="",xlab="",ylab="",
#                    label = "none") + # hide individual labels 
#   theme_cowplot() + 
#   theme(axis.text.y   = element_blank(),
#         axis.text.x   = element_blank(),
#         panel.border = element_rect(colour = "black", fill=NA, size=1))

pPCA_C <- plot_grid(p1,labels="A")  
pPCA_C

#### Supplemental Figures to PCA ####
# Function for generating three panel contribution plot for the different PCs
suppContriFig <- function(x){
  temp <- list()
  for(i in 1:3){
    temp[[i]] <- fviz_contrib(x, choice = "var", axes = i, top = 10,color = "black",fill = "grey",
                                          sort.val = "none")
  }
  return(temp)
}
# Control
contribute_C_fig <- suppContriFig(pca_C_obj) 
plot_grid(plotlist =contribute_C_fig, labels = c('A', 'B','C'), label_size = 12,nrow = 3)
p1_amongTrt
# Among Treatments
contribute_amongTrt_fig <- suppContriFig(p1_amongTrt) 
plot_grid(plotlist =contribute_amongTrt_fig, labels = c('A', 'B','C'), label_size = 12,nrow = 3)
#### Correlation plots between gene expression and DNA methylation ####
#Function for creating density coloring on scatterplots
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

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
plot_grid(plotlist=pL_mean,labels = namePlot)

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
plot_grid(plotlist=pL_cv,labels = namePlot)
# Development 
right <- plot_grid(pL_mean[[2]],pL_cv[[2]],labels=c("B","C"),ncol=1)
plot_grid(pPCA_C,right)


### Among Treatments ###
# CV vs diff methylation amount treatment
pL_Trtcv_9 <- list()
for(i in 1:length(dats)){
  df <- get(dats[1])
  df$density <- get_density(c(log(df$gene_cv_9E)-log(df$gene_cv_9C)),c(df$Mean_9E-df$Mean_9C)*100, n = 100)
  temp <- ggplot(df,aes(x=c(log(gene_cv_9E)-log(gene_cv_9C)),y=c(df$Mean_9E-df$Mean_9C)*100)) + 
    geom_point(aes(color = density)) + scale_color_viridis() +
    geom_smooth(method=lm,colour="orange",size=2) +
    labs(x="Delta log(CV_trt_tp9)",y="Delta DNA Methylation D9 (%)") +
    theme_cowplot()
  pL_Trtcv_9[[i]]<- temp
}


pL_Trtcv_80 <- list()
for(i in 1:length(dats)){
  df <- get(dats[2])
  temp <- rowMeans(cbind(df$Mean_80E,df$Mean_80C))
  df$diff_std <- c(df$Mean_80E-df$Mean_80C)/temp
  df$density <- get_density(c(log(df$gene_cv_80E)-log(df$gene_cv_80C)),
                            df$diff_std*100, n = 100)
                            #c(df$Mean_80E-df$Mean_80C)*100, n = 100)
  temp <- ggplot(df,aes(x=c(log(gene_cv_80E)-log(gene_cv_80C)),y=diff_std*100)) + 
    geom_point(aes(color = density)) + scale_color_viridis() +
    geom_smooth(method=lm,colour="orange",size=2) +
    labs(x="Delta log(CV_trt_tp80)",y="Delta DNA Methylation D80 (%)") +
    theme_cowplot()
  pL_Trtcv_80[[i]]<- temp
}

plot_grid(pL_Trtcv_9[[2]],pL_Trtcv_80[[2]])

  ## Model
  # Mean Methylation vs. Mean Gene Expression
  summary(lm(mat_20$gene_mean~mat_20$Mean_9C))
  # Mean Methylation vs. CV Gene Expression
  summary(lm(mat_20$gene_cv_9C~mat_20$Mean_9C))

#### DAPC DNAm vs Gene Expression ####

plotDAPC <- function(x,y,metaD){
  x_dapc <- dapc(t(x),metaD$Treatment,n.pca=7,n.da=1)
  y_dapc <- dapc(t(y),metaD$Treatment,n.pca=7,n.da=1)
  comb <- data.frame(Treatment=metaD$Treatment,
                     DNAm=y_dapc$ind.coord[,1],
                     Gene=x_dapc$ind.coord[,1])
  out <- ggplot(comb,aes(x=Gene,y=DNAm)) + geom_point() + 
          theme_cowplot() + 
          geom_smooth(method=lm) + 
          labs(x="Gene Body DAPC Coordinate",y="DNA Methylation DAPC Coordinate",title=paste0("Day ",unique(metaD$Day)[1])) +
          theme(plot.title = element_text(hjust = 0.5))
  return(out)
} 
# Two days to loop through
pL_dapcComparison <- list()
for(i in 1:2){
  meta_temp <- meta[meta$Day == unique(meta$Day)[i],]
  ge_counts_temp <- ge_counts[,meta$Day == unique(meta$Day)[i]]
  dnam_temp <- dnam[,meta$Day == unique(meta$Day)[i]]
  pL_dapcComparison[[i]] <- plotDAPC(ge_counts_temp,dnam_temp,meta_temp)
}
plot_grid(plotlist = pL_dapcComparison,ncol=1)

### Multi Plot ####
right_top <- plot_grid(pL_mean[[2]],pL_cv[[2]],labels=c("B","C"),ncol=1)
top <-  plot_grid(pPCA_C,right_top)
top
ggsave(top,filename="results/figures/Figure5_partial_topDevelopmental.png")
ggsave(top,filename="results/figures/Figure5_partial_topDevelopmental.pdf")
bottom <- plot_grid(pF,right_top,right_top,ncol=3,rel_widths = c(2,1,1))
final <- plot_grid(top,bottom,nrow=2)
final

#### Extra ####

## Difference in Methyl and Expression Among Treatments
plot(c(mat_20$Mean_9E-mat_20$Mean_9C)*100~c(mat_20$gene_mean_9E-mat_20$gene_mean_9C),
     xlab= "logfold GE Change at tp9", ylab="Diff. Methylation at tp9")
plot(c(mat_20$Mean_9E-mat_20$Mean_9C)~c(mat_20$gene_cv_mean_AmongTrt),xlim=c(-1,1))

mat_order <- mat[order(mat$diff_tp9_Trt),]
plot(mat_order$diff_tp9_Trt*100~c(1:length(mat$diff_tp9_Trt)))
matTP9 <-mat_order 
matTP9$pos <- c(1:length(mat$diff_tp9_Trt))
mat_sig_tp9 <- matTP9[!is.na(matTP9$DML_tp9_count),]

points(mat_sig_tp9$diff_tp9_Trt[order(mat_sig_tp9$diff_tp9_Trt)]*100~mat_sig_tp9$pos,col="red",pch=16)

plot(mat$diff_tp9_Trt~mat$cov5_count)
