#### Differential Gene Expresison ####

## packages
library(matrixStats,quietly = TRUE)
library(edgeR,quietly = TRUE)
library(limma,quietly = TRUE)
library(fdrtool,quietly = TRUE)
library(dplyr)
library(ggplot2)

#### Data ####
wd <- "~/Github/AE17_Cvirginica_MolecularResponse"
setwd(wd)
# This should be set to the path for the local version of the `2017OAExp_Oysters` github repo.
ge <- readRDS("data/Analysis/RNA_gene_postVoomAndNormalization_DGEListObj.RData")
#ge <- readRDS("data/Analysis/RNA_gene_postVoomAndNormalization_DGEListObj_remove17005.RData")
bioGenes <- read.csv("data/RNAseq/Target_BiomineralizationGenes.csv",stringsAsFactors = FALSE)
bioGenesV2 <- read.csv("data/RNAseq/Target_BiomineralizationGenes_V2.csv",stringsAsFactors = FALSE)
## Meta ##
meta <- readRDS("data/meta/metadata_20190811.RData")
meta$sampl_nameSimple <- substr(meta$sample_name,start = 4,stop=9)
#Create new factor levels (one for each level combination)
meta$SFVrn <- as.factor(paste0("D",meta$SFV))
meta$Sample_Index <- as.factor(meta$sample_index)
meta$TankID <- as.factor(meta$tankID)

# Remove 17005 since it appears as an outlier on PCA
# This was also run with 17005 and it doesn't impact interpretation
ge <- ge[,colnames(ge) != "17005"]
meta <- meta[meta$ID != "17005",]

# Performed differential expression with just biomineralization genes to confirm
# multiple hyp. corrections didn't impact the potential for significant genes after
# p adjustment.
# Result: It ultimately didn't alter the outcome.
#ge2 <- ge[which(row.names(ge) %in% bioGenes$Location),]
#ge <- ge2

## Design Matrix ##
design <- model.matrix(~0+SFVrn,data=meta) # 0+ is needed here otherwise the first level defaults to 1.
#Rename columns
colnames(design) <- levels(meta$SFVrn)

## Contrast Matrix ##
contr_mat <- makeContrasts(
  CvE_D9 = D09.2800-D09.400,
  CvE_D80 = D80.2800-D80.400,
  C_D9vD80 = D09.400-D80.400,
  Time = ((D09.2800-D09.400)- (D80.2800-D80.400))/2,
  Treatment = ((D09.2800+D80.2800)-(D09.400+D80.400))/2,
  levels=design
)

#### Analysis ####

#### Identify correlation between factors in design contrasts with blocking factor ####
ge_corr <- duplicateCorrelation(ge, design, block = meta$tankID)

#### Fitting Model ####
lmf_ge_corr <- lmFit(ge, design,
                        block = meta$tankID,
                        correlation = ge_corr$consensus.correlation)
lmf_ge_corr_V2 <- lmFit(ge, design)
# Ran the fit with and without correlation correction and outcome is same.

#Refitting Option 1 with contrasts
ge_contr <- contrasts.fit(lmf_ge_corr,contr_mat)
ge_contr_V2 <- contrasts.fit(lmf_ge_corr_V2,contr_mat)
##Run empiricial bayes protocol
ge_bayes <- eBayes(ge_contr,robust=TRUE)
ge_bayes_V2 <- eBayes(ge_contr,robust=TRUE)

## Top Candidates ##
## Criteria ( adj p value < 0.05 and lfc > 2)

top.table <- topTable(ge_bayes,number = Inf) 
head(top.table,5)

# Day 9 treatment comparison
top.table_d9 <- topTable(ge_bayes,lfc = 2,number = Inf,coef = 1)
top.table_V2_d9 <- topTable(ge_bayes_V2,lfc = 2,number = Inf,coef = 1) 
head(top.table_d9,5)
head(top.table_V2_d9,5)

# Day 80 treatment comparison
top.table_d80 <- topTable(ge_bayes,lfc = 2,number = Inf,coef = 2)
head(top.table_d80,5)
# No differentially expressed genes!

## Volcano Plots ##
coef_plot <- c(1,2) # Comparing treatments for each time point
par(mfrow=c(length(coef_plot),1))
for(i in coef_plot)volcanoplot(ge_bayes,coef = i,main=paste("Contrast",colnames(ge_bayes$contrasts)[i]))
par(mfrow=c(1,1))

## Manual p value calculation for each contrast using FDRtools and BH adjustment ##
# Confirmation using alternative multi-hyp correction approach (FDRtools and the t statistic)
colnames(ge_bayes$t)
min(fdrtool(ge_bayes$t[,1])$lfdr) #CvE_D9
min(fdrtool(ge_bayes$t[,2])$lfdr) #CvE_D80
min(fdrtool(ge_bayes$t[,3])$lfdr) #C_D9vD80
min(fdrtool(ge_bayes$t[,4])$lfdr) #Time
min(fdrtool(ge_bayes$t[,5])$lfdr) #Treatment
# Slight different minimum p values, but still nothing significant!

# Summary of all genes
top.table$Location <- row.names(top.table)
colSel <- c("Location","AveExpr","CvE_D9","CvE_D80","C_D9vD80",
            "Time","Treatment","adj.P.Val","predict")
tS_sub <- subset(top.table,select=colSel)
#write.csv(tS_sub,"paste0(wd,"/"data/Analysis/gene_diffExpression_all.csv")

## Summary of biomineralization genes

tS_bioGenes <- left_join(bioGenes,top.table,by=c("Location"))
tS_bioGenes <- tS_bioGenes[!is.na(tS_bioGenes$AveExpr),]
colSel <- c("Location","AveExpr","CvE_D9.y","CvE_D80.y","C_D9vD80.y",
            "Time","Treatment","adj.P.Val","predict")
tS_bioGenes_sub <- subset(tS_bioGenes,select=colSel)


tS_CA_sub <- subset(tS_CA,select=colSel)





library(reshape2)
library(cowplot)
library(matrixStats)
library("RColorBrewer")
pal <- brewer.pal(n = 12, name = 'Paired')
col_perm <- c(pal[1:2],pal[5:6],pal[12])

ge_df <- data.frame(Location=rownames(ge$E),ge$E,stringsAsFactors = FALSE)
colnames(ge_df)[-1] <- c(colnames(ge$E))

meta_ID <- data.frame(ID=meta$ID,SFV=meta$SFV,Treatment=meta$Treatment,Time=meta$Time)

# Original
x <- tS_bioGenes
# Version 2
x <- bioGenesV2
# Target Genes primarily associated with biomineralization
target_genes <- c("CA","CATP","NHE","NCX","BT")
x2 <- bioGenesV2[bioGenesV2$CV_GeneCode %in% target_genes,]

boxSummary_list <- createBoxplotSummary(ge_df,x,meta)
plot_grid(plotlist = boxSummary_list,ncol = 3)
  plot_grid(plotlist=boxSummary_list[c(1:12)],ncol = 3)
plot_grid(plotlist=boxSummary_list[c(13:length(boxSummary_list))],ncol = 3)

box_list <- createBoxplot(ge_df,x2,meta)
plot_grid(box_list[[2]],ncol = 1)
plot_grid(plotlist = box_list,ncol = 1)

createBoxplotSummary <-  function(ge_df=ge_df,x=yS_bioGenes,meta=meta){
  x <- x[!is.na(x$CV_GeneCode),]
  plot_list <- list()
  counter <- 1
  for(i in 1:length(unique(x$CV_GeneCode))){
    y <- x[x$CV_GeneCode == unique(x$CV_GeneCode)[i],]
    y <- y[!is.na(y$Location),]
    yy <- left_join(y,ge_df,by="Location")
    if(sum(!is.na(yy$`17007`))>0){
      yy <- yy[!is.na(yy$`17007`),]
      out <- as.data.frame(t(yy[,-c(1:8)]))
      out$ID <- rownames(out)
      if(ncol(out)> 2){
        out_summary <- data.frame(ID=out$ID,value=rowMeans(out[,-ncol(out)]))
      }else{
        out_summary <- data.frame(ID=out$ID,value=out[,1])
      }
      
      out_all_summary <- left_join(meta_ID,out_summary,by="ID")
      
      temp <- ggplot(out_all_summary,aes(x=SFV,
                                         y=value,
                                         fill=SFV)) + 
        geom_boxplot() + 
        scale_fill_manual(values = col_perm[c(1,2,3,4)])+
        theme_cowplot() + 
        labs(y= "GE (log2 - cpm)",
             x = "",
             title = paste0(y$CV_GeneCode," (n = ",nrow(yy),")")) +
        theme(plot.title = element_text(hjust=0.5),
              axis.title.x = element_blank(),
              legend.position = "none",
              axis.text.x = element_blank()) +
        scale_x_discrete(limits=c("09.400","09.2800","80.400","80.2800"),
                         labels=c("Day 9\nControl",
                                  "Day 9\nHigh OA",
                                  "Day 80\nControl",
                                  "Day 80\nHigh OA"))
      plot_list[[counter]] <- temp
      counter<-counter+1
    }
  }
  temp
  return(plot_list)
}
createBoxplot <- function(ge_df=ge_df,x=yS_bioGenes,meta=meta){
  plot_list <- list()
  for(i in 1:length(unique(x$CV_GeneCode))){
    y <- x[x$CV_GeneCode == unique(x$CV_GeneCode)[i],]
    y <- y[!is.na(y$Location),]
    yy <- left_join(y,ge_df,by="Location")
    out <- as.data.frame(t(yy[,-c(1:8)]))
    out$ID <- rownames(out)
    meta_ID_simple <- data.frame(ID=meta$ID,SFV=meta$SFV)
    out_all <- left_join(meta_ID_simple,out)
    out_long <- melt(out_all, id.var = c("SFV","ID"))
    temp <- ggplot(out_long,aes(x=interaction(SFV),y=value,fill=SFV)) + 
      facet_grid(.~variable) + 
      geom_boxplot() + 
      scale_fill_manual(values = col_perm[c(1,2,3,4)])+
      theme_cowplot() +
      labs(y= "GE (log2 - cpm)",
           x = "",
           title = paste0(y$CV_GeneCode," (n = ",nrow(y),")")) +
      theme(plot.title = element_text(hjust=0.5),
            axis.title.x = element_blank(),
            legend.position = "none",
            axis.text.x = element_blank()) +
      ylim(-5,12) +
      scale_x_discrete(limits=c("09.400","09.2800","80.400","80.2800"),
                       labels=c("Day 9\nControl",
                                "Day 9\nHigh OA",
                                "Day 80\nControl",
                                "Day 80\nHigh OA"))
    plot_list[[i]] <- temp
  }
return(plot_list)
}

#write.csv <- (tS_bioGenes_sub,"paste0(wd,"/"data/Analysis/gene_diffExpression_biomineralizationGenes.csv")

#  Notes: 
# $coefficients in the lmFit object are the log2 expression for each treatment x time level
# $coefficients in the eBayes objects are the log2 fold change (calculated at the difference, A-B) for each
#               contrast

#### Additional analyses not in the manuscript ####

#### Alternative thresholds for significant ##

## We have a number of genes that exhibit some pretty large (> abs(2) fold changes) and relatively
#  large -log10(pvalues), but these were not significant after correction.

## Lets take a look at some with larger fold changes (>abs(2)) and -log10(p.values) > 2.5
ge_bayes_top <- ge_bayes[abs(ge_bayes$coefficients[,1]) > 2 & -log10(ge_bayes$p.value[,1]) > 2.5,] 
ge_bayes_top$coefficients
# LOC111137004 : bridging integrator 2-like - downregulated in D9_2800
# LOC111122726 : heatshock protein - upregulated in D9_2800
# LOC111121227 : methylthioribose kinase-like - upregulated in 2800 treatment
# Visualization of top candidates
target <- ge$E[row.names(ge$E) == "LOC111137004",]
target.df <- data.frame(ge=target,trt=meta$SFVrn)
ggplot(target.df,aes(y=target,x=trt)) + geom_boxplot()
# Looks like an interesting differential expression candidate but not ultimately significant.

#### Correlating Gene Expression with EPF pH ##
  ## Here we use the simple linear model from the limma package 
  ## to explore the relationship between expression and EPF pH (or delta EPF pH)
 
## Extra pallial fluid pH
design_PH <- model.matrix(~epf_pH,data=meta)
design_PH <- model.matrix(~diff_pH,data=meta)
ge_corr_PH <- duplicateCorrelation(ge, design_PH, block = meta$tankID)
lmf_ge_corr_PH <- lmFit(ge, design_PH,
                         block = meta$tankID,
                         correlation = ge_corr_PH$consensus.correlation)
contr_PH <- contrasts.fit(lmf_ge_corr_PH, coef = 2) # test "pH" coefficient
bayes_PH <- eBayes(contr_PH)
top.table.pH <- topTable(bayes_PH, sort.by = "P", n = Inf)
head(top.table.pH, 20)
## delta EPF pH (EPF pH - Environmental pH)
design_dPH <- model.matrix(~diff_pH,data=meta)
ge_corr_dPH <- duplicateCorrelation(ge, design_dPH, block = meta$tankID)
lmf_ge_corr_dPH <- lmFit(ge, design_dPH,
                     block = meta$tankID,
                     correlation = ge_corr_dPH$consensus.correlation)
contr_dPH <- contrasts.fit(lmf_ge_corr_dPH, coef = 2) # test "pH" coefficient
bayes_dPH <- eBayes(contr_dPH)
top.table.dPH <- topTable(bayes_dPH, sort.by = "P", n = Inf)
head(top.table, 20)

# Nothing is significant!
