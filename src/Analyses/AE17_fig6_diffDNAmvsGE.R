#### Correlated DNA methylation and gene expression responses to OA ####
# Figures from this code : 6
  
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
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
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
ge<-readRDS("results/RNA/RNA_gene_postVoomAndNormalization_DGEListObj.RData")
ge_counts <- ge$E # Extract just the counts
ge_counts <- ge_counts[,colnames(ge_counts) != "17099"] # Remove problematic individual 
ge_counts <- ge_counts[,colnames(ge_counts) != "17005"] # Remove problematic individual 

# Read in DNA Methylation Data for Gene Bodys
dnam <- readRDS("results/DNAm/DNAm_20200202_AllCountsList_cov5_byFeature.RData")
dnam <-dnam$beta$gene # Only going to look at the beta (dna methylation value) for genes
dnam <- as.matrix(dnam[,2:ncol(dnam)]) # Change data into matrix 
dnam[is.na(dnam)] <- 0
class(dnam) <- "numeric"
dnam <- dnam[,colnames(dnam) != "17005"] # Remove problematic individual 

# Gene Atribute, Expression and Methylation Summary Table
mat <- readRDS("results/Multi/Multi_geneSummaryReduced.RData")
# We have overlapping expression and DNA Methylation data for 9626 genes

# Read in list of significant DMLs
dmls <- read.csv("results/DNAm/DNAm_methylkit_dm50_q01_sigSummaryTable.csv")
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

dmls_9_join$diff.cv <- c(dmls_9_join$gene_cv_9E-dmls_9_join$gene_cv_9C)
dmls_80_join$diff.cv <- c(dmls_80_join$gene_cv_80E-dmls_80_join$gene_cv_80C)

dmls_9_join$methylDirection <- ifelse(dmls_9_join$meth.diff < 0,
                                      "Hypomethylation",
                                      "Hypermethylation")
dmls_80_join$methylDirection <- ifelse(dmls_80_join$meth.diff < 0,
                                       "Hypomethylation",
                                       "Hypermethylation")

#### Gene Expression vs. DML change ####
p1<-ggplot(dmls_9_join,
           aes(x=meth.diff,y=gene_diff_tp9_Trt,colour=methylDirection)) + 
  geom_smooth(method=lm,colour="orange",size=2,na.rm = TRUE) +
  geom_point(size=3) + xlim(-100,100) + ylim(-2,2) + #ylim(-1.1,1.1) + 
  theme_cowplot() + scale_color_manual(values=pal[c(3,1)]) +
  labs(y=bquote("Gene expression (FoldChange"~Log[2]~")"),x="Methylation change in DMLs (difference %)",title="Day 9") +
  theme(plot.title = element_text(hjust=0.5))
p1_nolab <- p1 + theme(axis.title = element_blank(),
                       legend.position = "none")
p1_nolab
p2<-ggplot(dmls_80_join[!is.na(dmls_80_join$gene_diff_tp80_Trt),],
           aes(x=meth.diff,y=gene_diff_tp80_Trt,colour=methylDirection)) + 
  geom_smooth(method=lm,colour="orange",size=2,na.rm = TRUE) +
  geom_point(size=3) + xlim(-100,100) + ylim(-2,2) + #ylim(-1.1,1.1) + 
  theme_cowplot() + scale_color_manual(values=pal[c(3,1)]) +
  labs(y=bquote("Gene expression (FoldChange"~Log[2]~")"),x="Methylation change in DMLs (difference %)",title="Day 80") +
  theme(plot.title = element_text(hjust=0.5))
p2_nolab <- p2 + theme(axis.title = element_blank(),
                       legend.position = "none")
p2_nolab

## Statistical Model ##
# Day 09 
#plot(lm(dmls_9_join$diff_tp9_Trt~dmls_9_join$meth.diff))
summary(lm(dmls_9_join$diff_tp9_Trt~dmls_9_join$meth.diff))
# Slope 0.0004903 
# P     0.000484
# R2    0.4061

# Day 80  
dml_80s <- dmls_80_join[!is.na(dmls_80_join$gene_diff_tp80_Trt),]
#plot(lm(dml_80s$diff_tp80_Trt~dml_80s$meth.diff))
# slight outlier but doesn't impact outcome all that much
summary(lm(dml_80s$diff_tp80_Trt~dml_80s$meth.diff))
# Slope 3.973e-04 
# P     3.77e-05
# R2    0.2762

#### Gene level summary difference ####
mat_final <- mat[mat$cov5_count/mat$all_count >= 0.2,]
day9_dnam_diff <- (mat_final$Mean_9E-mat_final$Mean_9C)*100
day80_dnam_diff <- (mat_final$Mean_80E-mat_final$Mean_80C)*100
day9_ge_diff <- mat_final$gene_diff_tp9_Trt
day80_ge_diff <- mat_final$gene_diff_tp80_Trt

gene_summary <- data.frame(d9_dnam=day9_dnam_diff,
                           d9_ge=day9_ge_diff,
                           d80_dnam=day80_dnam_diff,
                           d80_ge=day80_ge_diff)

gene_summary$D9_density <- get_density(gene_summary$d9_ge, gene_summary$d9_dnam, n = 100)
gene_summary$D80_density <- get_density(gene_summary$d80_ge, gene_summary$d80_dnam, n = 100)

## Day 9 - Gene level summary difference
pG_9<-ggplot(gene_summary,
           aes(y=day9_ge_diff,x=d9_dnam)) + 
  geom_point(size=3,aes(color=D9_density)) + 
  xlim(-25,25) + ylim(-2,2) +
  scale_color_viridis() +
  geom_smooth(method=lm,colour="orange",size=2,na.rm = TRUE) +
  theme_cowplot() + 
  labs(y=bquote("Gene expression (FoldChange"~Log[2]~")"),x="Mean gene methylation (difference %)",title="Day 9") +
  theme(plot.title = element_text(hjust=0.5))
pG_9
pG_9_nolab <- pG_9 + theme(axis.title = element_blank(),
                       legend.position = "none")
pG_9_nolab

## Day 80 - Gene level summary difference
pG_80<-ggplot(gene_summary,
           aes(y=day80_ge_diff,x=d80_dnam)) + 
  geom_point(size=3,aes(color=D80_density)) + 
  xlim(-25,25) + ylim(-2,2) +
  scale_color_viridis() +
  geom_smooth(method=lm,colour="orange",size=2,na.rm = TRUE) +
  theme_cowplot() + 
  labs(y=bquote("Gene expression (FoldChange"~Log[2]~")"),x="Mean gene methylation (difference %)",title="Day 80") +
  theme(plot.title = element_text(hjust=0.5))
pG_80
pG_80_nolab <- pG_80 + theme(axis.title = element_blank(),
                       legend.position = "none")
pG_80_nolab

plot_grid(pG_9,pG_80,nrow=1)

## Statitical models ###
#Day 9
summary(lm(gene_summary$d9_ge~gene_summary$d9_dnam))
# slope 0.012131
# R2  0.01243
# P 1.15e-11

# Day 80
summary(lm(gene_summary$d80_ge~gene_summary$d80_dnam))
# Slop 0.0137262
# P   3.1e-13
# R2  0.01437

#### Final Figure ####
plot1 <- plot_grid(p1_nolab,p2_nolab,ncol=2)
plot2 <- plot_grid(pG_9_nolab,pG_80_nolab,ncol=2)
x.grob1 <- textGrob("Differential Methylated Loci (difference %)", 
                    gp=gpar(col="black", fontsize=15))
x.grob2 <- textGrob("DNA Methylation (difference %)", 
                    gp=gpar(col="black", fontsize=15))
y.grob1 <- textGrob("Gene Expression (log2 fold change)", 
                    gp=gpar( col="black", fontsize=15), rot=90)

save1<-grid.arrange(arrangeGrob(plot1,bottom = x.grob1))
save2<-grid.arrange(arrangeGrob(plot2, bottom = x.grob2))
plot3<-plot_grid(save2,save1,nrow=2,labels=c("A","B"))
saveG<-grid.arrange(arrangeGrob(plot3, left = y.grob1))

