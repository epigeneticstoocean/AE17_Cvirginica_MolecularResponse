#### Differential Gene Expresison ####

## packages
library(matrixStats,quietly = TRUE)
library(edgeR,quietly = TRUE)
library(limma,quietly = TRUE)
library(fdrtool,quietly = TRUE)

#### Data ####
wd <- "~/Github/AE17_Cvirginica_MolecularResponse"
# This should be set to the path for the local version of the `2017OAExp_Oysters` github repo.
ge <- readRDS(paste0(wd,"/data/Analysis/RNA_gene_postVoomAndNormalization_DGEListObj.RData"))

## Meta ##
meta <- readRDS(paste0(wd,"/data/meta/metadata_20190811.RData"))
meta$sampl_nameSimple <- substr(meta$sample_name,start = 4,stop=9)
#Create new factor levels (one for each level combination)
meta$SFVrn <- as.factor(paste0("D",meta$SFV))
meta$Sample_Index <- as.factor(meta$sample_index)
meta$TankID <- as.factor(meta$tankID)

## Design Matrix ##
design <- model.matrix(~0+SFVrn,data=meta) # 0+ is needed here otherwise the first level defaults to 1.
#Rename columns
colnames(design) <- levels(meta$SFVrn)

## Contrast Matrix ##
contr_mat <- makeContrasts(
  CvE_D9 = D09.2800-D09.400,
  CvE_D80 = D80.2800-D80.400,
  C_D9vD80 = D09.400-D80.400,
  Time = (D09.2800-D09.400)- (D80.2800-D80.400),
  Treatment = (D09.2800-D80.2800)- (D09.400-D80.400),
  levels=design
)

contr_mat <- makeContrasts(
  CvE_D9 = D09.2800-D09.400,
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
head(ge_bayes$coefficients)
#  Notes: 
# $coefficients in the lmFit object are the log2 expression for each treatment x time level
# $coefficients in the eBayes objects are the log2 fold change (calculated at the difference, A-B) for each
#               contrast

## Volcano plots for visualizing differential gene expression
volcanoplot(ge_bayes)
# Equivalent code in ggplot
temp <- data.frame(FoldChange=as.numeric(ge_bayes$coefficients),p.value=-log10(as.numeric(ge_bayes$p.value)))
ggplot(temp,aes(x=FoldChange,y=p.value)) + geom_point()
## We have a number of genes that exhibit some pretty large (> abs(2) fold changes) and relatively
 # large -log10(pvalues), but these were not significant after correction.

## Lets take a look at some with larger fold changes (>abs(2)) and -log10(p.values) > 2.5
ge_bayes_top <- ge_bayes[abs(ge_bayes$coefficients) > 2 & -log10(ge_bayes$p.value) > 2.5,] 
ge_bayes_top$coefficients
# LOC111137004 : bridging integrator 2-like - downregulated in D9_2800
# LOC111122726 : heatshock protein - upregulated in D9_2800
# LOC111121227 : methylthioribose kinase-like - upregulated in 2800 treatment
target <- ge$E[row.names(ge$E) == "LOC111121227",]
target.df <- data.frame(ge=target,trt=meta$SFVrn)
ggplot(target.df,aes(y=target,x=trt)) + geom_boxplot()

## Saving Files
#saveRDS(bayes_gene_a2_o1_contr,paste0(wd,"/gene_EBayesObj.RData"))
