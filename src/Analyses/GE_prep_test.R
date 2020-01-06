#### Differential Expression Analyses #####

## Packages 
library(matrixStats,quietly = TRUE)
library(edgeR,quietly = TRUE)
library(dplyr,quietly = TRUE)
library(FSA,quietly = TRUE) # to use w arguement in hist()
library(ggplot2)
library(cowplot)
library(limma)

#### Data #### 
wd <- "/home/downeyam/Github/2017OAExp_Oysters"
# This should be set to the path for the local version of the `2017OAExp_Oysters` github repo.

## Meta Data for the Model
#model<-readRDS("/home/downeyam/Github/2017OAExp_Oysters/input_files/meta/metadata_20190811.RData")
# Gene
gc_diff <- readRDS(paste0(wd,"/results/Transcriptomic/gene_EBayesObj.RData"))
# Transcripts
tc_diff <- readRDS(paste0(wd,"/results/Transcriptomic/transcript_EBayesObj.RData"))
# Biomineralization Genes
tl <- read.csv(paste0(wd,"/input_files/RNA/references/Target_BiomineralizationGenes.csv"))
## Top table of Genes
g_top <- topTable(gc_diff,number = Inf)
## Top table of Transcripts
t_top <- topTable(tc_diff,number=Inf)

### CvE on Day 9 : CvE on Day 80  
```{r}
par(mfrow=c(1,1))
tl$Location <- as.character(tl$Location)
gc_diff_bio <- gc_diff[which(!is.na(match(gc_diff$genes$GENEID,tl$Location))),]

plot(gc_diff$coefficients[,1]~gc_diff$coefficients[,2],
     col=adjustcolor("red",alpha.f = 0.4),
     xlim=c(-8,8),ylim=c(-8,8),
     xlab="Log2 Fold Change : Day 80",
     ylab="Log2 Fold Change : Day 9")
abline(h=0,col="black")
abline(v=0,col="black")
points(gc_diff_bio$coefficients[,1]~gc_diff_bio$coefficients[,2],
       col=adjustcolor("blue",alpha.f=0.7),pch=16)

dat <- readRDS("/home/downeyam/Github/2017OAExp_Oysters/results/Transcriptomic/EBayesObj_gene_withIndWeights_filterApproach2_plannedContrastMatrix.RData")
cds <- readRDS(paste0(wd,"/input_files/RNA/references/CDS_wGoTerms_GeneLOC.RData"))
# Extract GO slim terms
go_tab_slim <- cds[,c(12,14)]
go_tab <- cds[,c(12,17)]
go_tab_unique <- go_tab[!duplicated(go_tab$gene_id),]
go_tab_unique$InterproScan.GO.ID <- as.character(go_tab_unique$InterproScan.GO.ID)
go_tab_unique$InterproScan.GO.ID[is.na(go_tab_unique$InterproScan.GO.ID)] <- "unknown"
go_tab_unique$InterproScan.GO.ID[go_tab_unique$InterproScan.GO.ID==""] <- "unknown"
#creat two datasets, one for each timepoint comparison
t9 <- data.frame(gene_id=as.character(rownames(dat$coefficients)),logfold=dat$coefficients[,1],stringsAsFactors = FALSE)
t80 <- data.frame(gene_id=as.character(rownames(dat$coefficients)),logfold=dat$coefficients[,2],stringsAsFactors = FALSE)
# Join cpg annotations with go slim terms
comb_9 <- left_join(t9,go_tab_unique,"gene_id")
comb_80 <- left_join(t80,go_tab_unique,"gene_id")
comb_9_red <- comb_9[c(1,3)]
comb_80_red <- comb_80[c(1,3)]
head(go_tab_unique)
write.table(x = go_tab_unique,file = "/home/downeyam/Github/2017OAExp_Oysters/results/Transcriptomic/GOMWU_LogFoldChangeGE_Goterms.tab",sep = "\t",row.names = FALSE)
write.csv(x = t9,file = "/home/downeyam/Github/2017OAExp_Oysters/results/Transcriptomic/GOMWU_LogFoldChangeGE_T9.tab",row.names = FALSE)
write.csv(x = t80,file = "/home/downeyam/Github/2017OAExp_Oysters/results/Transcriptomic/GOMWU_LogFoldChangeGE_T80.tab",row.names = FALSE)
