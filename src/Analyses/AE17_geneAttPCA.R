#### Gene Expression, Methylation, and Gene Attribute PCA Analysis ####

## package
knitr::opts_chunk$set(echo = TRUE)
library(rtracklayer) # on bioconducter
library(dplyr)
library(edgeR)
library(factoextra)
library(cowplot)

#### Data ####
## Gene Expression
RSEM <-  readRDS("RSEM_gene_Summary.Rdata")
# Separate out RSEM counts and rename rows with LOC ID
rsem_c <- RSEM$Count # Stick with raw gene counts  (not TPM or some other normalizaed count)
rm(RSEM)
## Gene Annotation
# Transcript file
tran <- readRDS("references/STAR_gnomon_tximportGeneFile.RData")
# Gene File
gene <- tran[!duplicated(tran$GENEID),]
gene$gene_length <- gene$stop-gene$start
## Meta Data 
meta <- readRDS("meta/metadata_20190811.RData")
meta$sampl_nameSimple <- substr(meta$sample_name,start = 4,stop=9)
#Create new factor levels (one for each level combination)
meta$SFVrn <- as.factor(paste0("D",meta$SFV))
meta$Sample_Index <- as.factor(meta$sample_index)
meta$TankID <- as.factor(meta$tankID)
## Data Manipulation 
# Order genes from annotation list to match the order in the count matrix
gene_order <- gene[order(gene$gene_id),]
identical(rownames(rsem_c),gene_order$gene_id) # TRUE confirms the order matches
# Relabel the rows of the count matrix with LOC ID
rownames(rsem_c) <- gene_order$GENEID
geneC_all <- cpm(round(rsem_c))
rm(rsem_c)
# CpGs with coverage of 5
cpgMean <-  read.table("DNAm/mean_betamean_CVTrtMeans_PerFeature_filteredBetameanGT05.txt"),header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
#Keep only the gene rows
cpgMean_all_gene <- cpgMean_all[cpgMean_all$V3 == "gene",]
cpgMean_gene <- cpgMean[cpgMean$V3 == "gene",]
# Remove the full datasets since they are large
# rm(cpgMean)
# rm(cpgMean_all)

## Exon Count per Gene 
cpgMean_all_exon <- cpgMean_all[cpgMean_all$V3 == "exon",]
gene_type_byExon <- sub(".*gbkey=(.*?);.*","\\1",cpgMean_all_exon$V9,perl=TRUE)
# Selecting out exons within genes we care about  (no t or rRNA)
tf_mRNA <- gene_type_byExon == "mRNA" | gene_type_byExon == "ncRNA" | gene_type_byExon == "misc_RNA" |  gene_type_byExon == "exon"
#Pulling out geneID 
gene_ID_inExons <- sub(".*GeneID:(.[0-9]{1,8}).*","\\1",cpgMean_all_exon$V9,perl=TRUE)
gene_ID_inExons_filtered <- as.numeric(gene_ID_inExons[tf_mRNA])
gene_ID_inExons_filtered <- paste0("LOC",gene_ID_inExons_filtered)
cpgMean_all_exon$GENEID[tf_mRNA] <- gene_ID_inExons_filtered
cpgMean_all_exon$GENEID[!tf_mRNA] <-"NA"
# Generate Exon counts
exon_counts <- data.frame(GENEID=names(table(cpgMean_all_exon$GENEID)),Exon_Count=as.vector(table(cpgMean_all_exon$GENEID)))
# Count number of CpGs in exons
exon_cpg_counts <- aggregate(V11~GENEID,FUN=sum,data=cpgMean_all_exon)
# Take mean MEthylation for cpgs in exons
cpgMean_all_exon$V10 <- as.numeric(cpgMean_all_exon$V10)
exon_cpg_means <- aggregate(V10~GENEID,FUN=mean,data=cpgMean_all_exon,na.rm=TRUE)
# Merge these piece of data together
exons <- left_join(exon_counts,exon_cpg_counts)
exons <- left_join(exons,exon_cpg_means)
#remove the na
exons <- exons[!exons$GENEID == "NA",]
colnames(exons)[3:4] <- c("exon_cpg_counts","exon_cpg_means")
## Extract Gene ID
# This is to merge DNA methylation data with gene expression data
#Total CpG data
gene_ID <- sub(".*GeneID:(.*?);.*","\\1",cpgMean_all_gene$V9,perl=TRUE)
LOC_ID <- paste0("LOC",gene_ID)
cpgMean_all_gene_ID <- data.frame(GENEID=LOC_ID,cpgMean_all_gene)
# Filtered data
gene_ID <- sub(".*GeneID:(.*?);.*","\\1",cpgMean_gene$V9,perl=TRUE)
LOC_ID <- paste0("LOC",gene_ID)
cpgMean_gene_ID <- data.frame(GENEID=LOC_ID,cpgMean_gene)

plot(cpgMean_all_gene_ID$V10~cpgMean_gene_ID$V10,
     xlim=c(0,1),ylim=c(0,1),
     xlab="Filter CpGs",
     ylab="Unfiltered CpGs",
     main="Gene Level DNA Methylation Percent")
abline(b=1,a=0,col="red",lwd=2)

#Creating new data frame with counts
cpg_final <- data.frame(cpgMean_gene_ID[,c(1,11:13)],
                        cpgMean_total=cpgMean_all_gene_ID[,11],
                        cpgCount_total=cpgMean_all_gene_ID[,12])
colnames(cpg_final)[2:4] <- c("cpgMean_filtered","cpgCV_filtered","cpgCount_filtered")
# Merging gene expression, gene leve dnam, and exon level dnam
initial <- left_join(gene,cpg_final)
ge_dnam <- left_join(initial,exons)
ge_dnamV2 <- ge_dnam[!is.na(ge_dnam$cpgCount_total),]
out <- match(ge_dnamV2$GENEID,rownames(geneC_all))
# Reorder and subset gene expression matrix to match dnam
geneC_new <- geneC_all[match(ge_dnamV2$GENEID,rownames(geneC_all)),]
#Confirm that the dnam matrix and gene matrix match
identical(rownames(geneC_new),ge_dnamV2$GENEID)
## Filtering by genes with low or no expression
#Genes Filtering (same as for differential expression analysis)
# Breaking down expression coverage by treatment*time combination
#Day 9 Trt 2800
keep_D9.2800 <- rowSums(geneC_new[,meta$SFVrn=="D09.2800"]>=1) >= 5
sum(keep_D9.2800)
#Day 9 Trt 400
keep_D9.400 <- rowSums(geneC_new[,meta$SFVrn=="D09.400"]>=1)>= 5
sum(keep_D9.400)
#Day 80 Trt 2800
keep_D80.2800 <- rowSums(geneC_new[,meta$SFVrn=="D80.2800"]>=1) >= 5
sum(keep_D80.2800)#
#Day 80 Trt 400
keep_D80.400 <- rowSums(geneC_new[,meta$SFVrn=="D80.400"]>=1) >= 5
sum(keep_D80.400)

keep_gene_a2 <- rowSums(cbind(keep_D9.2800,keep_D9.400,
                              keep_D80.2800,keep_D80.400)) >= 1
# Filter low coverage genes (remember data is already in cpms)
geneC_final <- geneC_new[keep_gene_a2, ]
# Filter dnam data similarly
ge_dnamV3 <- ge_dnamV2[keep_gene_a2,]
## Summarize Gene Expression
cv <- function(x)(sd(x)/mean(x))
gene_Mean <-  rowMeans(geneC_final)
gene_cv <-  apply(geneC_final,1,cv)
gene_sum_all <- data.frame(gene_Mean,gene_cv)

gene_Mean_9C <- rowMeans(geneC_final[,meta$SFVrn=="D09.400"])
gene_cv_9C <-  apply(geneC_final[,meta$SFVrn=="D09.400"],1,cv)
gene_sum_9C <- data.frame(gene_Mean_9C,gene_cv_9C)

gene_Mean_9E <- rowMeans(geneC_final[,meta$SFVrn=="D09.2800"])
gene_cv_9E <-  apply(geneC_final[,meta$SFVrn=="D09.2800"],1,cv)

gene_Mean_80C <- rowMeans(geneC_final[,meta$SFVrn=="D80.400"])
gene_cv_80C <-  apply(geneC_final[,meta$SFVrn=="D80.400"],1,cv)

gene_Mean_80E <- rowMeans(geneC_final[,meta$SFVrn=="D80.2800"])
gene_cv_80E <-  apply(geneC_final[,meta$SFVrn=="D80.2800"],1,cv)

trt_cv <- apply(cbind(gene_Mean_9C,gene_Mean_9E,gene_Mean_80C,gene_Mean_80E),1,cv)
diff <- gene_Mean_9C-gene_Mean_9E

ge_sum <- data.frame(gene_Mean,trt_cv)
colnames(ge_sum) <- c("Expression","Expression_CV")
hist(ge_sum$Expression_CV)
## Format data for PCA
col <- c("gene_length","cpgMean_filtered","cpgCV_filtered","cpgCount_total","Exon_Count") 
cpg_pca <- subset(ge_dnamV3,select=col)
colnames(cpg_pca) <- c("gene_length","Methylation","Methylation_CV","CpGs","Exons")
cpg_pca$Methylation <- as.numeric(as.character(cpg_pca$Methylation))
cpg_pca$Methylation_CV<- as.numeric(as.character(cpg_pca$Methylation_CV))
## Futher reducing the dnam matrix and gene count matrix down to size.
final_pca <- data.frame(ge_sum,cpg_pca)
final_pca <- final_pca[!is.na(final_pca$Methylation | final_pca$Methylation_CV),]
# log2 transform
final_pca2 <- final_pca
for(i in 1:ncol(final_pca)){
  if(i >= 3 & i < 6){final_pca2[,i] <- asin(final_pca[,i])
  }else{final_pca2[,i] <- log2(final_pca[,i])}
}
final_pca3 <- scale(final_pca2)

#### Analysis ####

# We are using log2 transformed and scaled variables (except methylation is not log2 transformed)

# Perform pca
pca_obj <- prcomp(final_pca3)
# Scree plot 
fviz_eig(pca_obj,addlabels = TRUE)

#### Figure 4 ####
# Black and white
p1<-fviz_pca_var(pca_obj,col.circle = "white",
                 title="",
                 repel = TRUE # Avoid text overlapping
) + theme_classic() +
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=15))
p2 <- fviz_pca_ind(pca_obj,title="",xlab="",ylab="",
                   label = "none") + # hide individual labels 
  theme_classic() + 
  theme(axis.text.y   = element_blank(),
        axis.text.x   = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

pF <- ggdraw(p1) +
  draw_plot(p2,0.19,.09,.33,.4) 
pF
draw_plot_label(
  c("A", "B"),
  c(0.05, 0.12),
  c(1, 0.45),
  size = 20
)

#### Supplemental Figure ####
# Contributions of variables to PC1
p3 <- fviz_contrib(pca_obj, choice = "var", axes = 1, top = 10,color = "black",fill = "grey",
                   sort.val = "none")
# Contributions of variables to PC2
p4 <- fviz_contrib(pca_obj, choice = "var", axes = 2, top = 10,color = "black",fill = "grey",
                   sort.val = "none")
# Contributions of variables to PC3
p5 <- fviz_contrib(pca_obj, choice = "var", axes = 3, top = 10,color = "black",fill = "grey",
                   sort.val="none")
fviz_contrib(pca_obj, choice = "var", axes = 3,
             top = 10,color = "black",fill = "grey",addLabels=TRUE)
pca_obj$rotation
plot_grid(p3,p4,p5, labels = c('A', 'B','C'), label_size = 12,nrow = 3)