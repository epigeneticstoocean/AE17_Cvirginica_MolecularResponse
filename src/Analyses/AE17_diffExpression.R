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
# This should be set to the path for the local version of the `2017OAExp_Oysters` github repo.
ge <- readRDS(paste0(wd,"/data/Analysis/RNA_gene_postVoomAndNormalization_DGEListObj.RData"))

bioGenes <- read.csv("data/RNAseq/Target_BiomineralizationGenes.csv",stringsAsFactors = FALSE)

## Meta ##
meta <- readRDS(paste0(wd,"/data/meta/metadata_20190811.RData"))
meta$sampl_nameSimple <- substr(meta$sample_name,start = 4,stop=9)
#Create new factor levels (one for each level combination)
meta$SFVrn <- as.factor(paste0("D",meta$SFV))
meta$Sample_Index <- as.factor(meta$sample_index)
meta$TankID <- as.factor(meta$tankID)

# Remove 17005 since it appears as an outlier on PCA
# This was also run with 17005 and it doesn't impact interpretation
ge <- ge[,colnames(ge) != "17005"]
meta <- meta[meta$ID != "17005",]

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
#Refitting Option 1 with contrasts
ge_contr <- contrasts.fit(lmf_ge_corr,contr_mat)
##Run empiricial bayes protocol
ge_bayes <- eBayes(ge_contr,robust=TRUE)
# Top Candidates
topTable(ge_bayes)
# No differentially expressed genes!

# Confirmation using alternative multi-hyp correction approach (FDRtools and the t statistic)
colnames(ge_bayes$t)
min(fdrtool(ge_bayes$t[,1])$lfdr) #CvE_D9
min(fdrtool(ge_bayes$t[,2])$lfdr) #CvE_D80
min(fdrtool(ge_bayes$t[,3])$lfdr) #C_D9vD80
min(fdrtool(ge_bayes$t[,4])$lfdr) #Time
min(fdrtool(ge_bayes$t[,5])$lfdr) #Treatment
# Still nothing significant!

# Summary of all genes
tS <- topTable(ge_bayes,number = 100000) # the 100,000 here just ensures all genes included
tS$Location <- row.names(tS)
colSel <- c("Location","AveExpr","CvE_D9","CvE_D80","C_D9vD80",
            "Time","Treatment","adj.P.Val","predict")
tS_sub <- subset(tS,select=colSel)
write.csv(tS_sub,paste0(wd,"/data/Analysis/gene_diffExpression_all.csv"))

## Summary of biomineralization genes
tS_bioGenes <- left_join(bioGenes,tS,by=c("Location"))
tS_bioGenes <- tS_bioGenes[!is.na(tS_bioGenes$AveExpr),]
colSel <- c("Location","AveExpr","CvE_D9.y","CvE_D80.y","C_D9vD80.y",
            "Time","Treatment","adj.P.Val","predict")
tS_bioGenes_sub <- subset(tS_bioGenes,select=colSel)
write.csv(tS_bioGenes_sub,paste0(wd,"/data/Analysis/gene_diffExpression_biomineralizationGenes.csv"))

#  Notes: 
# $coefficients in the lmFit object are the log2 expression for each treatment x time level
# $coefficients in the eBayes objects are the log2 fold change (calculated at the difference, A-B) for each
#               contrast

## Volcano plots for visualizing differential gene expression
volcanoplot(ge_bayes)
# Equivalent code in ggplot
temp <- data.frame(FoldChange=as.numeric(ge_bayes$coefficients),p.value=-log10(as.numeric(ge_bayes$p.value)))
ggplot(temp,aes(x=FoldChange,y=p.value)) + geom_point()


#### Alterntive thresholds for significant (not used in manuscript) ####

## We have a number of genes that exhibit some pretty large (> abs(2) fold changes) and relatively
#  large -log10(pvalues), but these were not significant after correction.

## Lets take a look at some with larger fold changes (>abs(2)) and -log10(p.values) > 2.5
ge_bayes_top <- ge_bayes[abs(ge_bayes$coefficients[,1]) > 2 & -log10(ge_bayes$p.value[,1]) > 2.5,] 
ge_bayes_top$coefficients
# LOC111137004 : bridging integrator 2-like - downregulated in D9_2800
# LOC111122726 : heatshock protein - upregulated in D9_2800
# LOC111121227 : methylthioribose kinase-like - upregulated in 2800 treatment
row.names(ge_bayes$coefficients)
ge_bayes$coefficients[row.names(ge_bayes$coefficients) == "LOC111099458",]
target <- ge$E[row.names(ge$E) == "LOC111129447",]
# Average Expression
mean(target)
target.df <- data.frame(ge=target,trt=meta$SFVrn)
ggplot(target.df,aes(y=target,x=trt)) + geom_boxplot()

## Saving Files
#saveRDS(bayes_gene_a2_o1_contr,paste0(wd,"/gene_EBayesObj.RData"))
