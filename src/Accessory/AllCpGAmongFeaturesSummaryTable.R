## Script creates a summary table of the number of cytosines from the different coverage lists 
 #(i.e. all CpGs, 5x coverage, etc) in each feature (i.e. exons, introns, etc). 

#### Data ####
## path
inputDir <- "~/Github/AE17_Cvirginica_MolecularResponse/data/MBDBS_seq"
outputDir <- "~/Github/AE17_Cvirginica_MolecularResponse/data/Analysis"

## Specific Feature coverage
setwd(inputDir)
# All cytosines in the genome
all_exon <- read.delim("20200129_CpG_all_Exon.txt",sep="\t",header = FALSE)
all_intron <- read.delim("20200129_CpG_all_Intron.txt",sep="\t",header = FALSE)
all_gene <- read.delim("20200129_CpG_all_gene.txt",sep="\t",header = FALSE)
all_intergenic <- read.delim("20200129_CpG_all_Intergenic.txt",sep="\t",header = FALSE)
# All cytosines with 5x coverage
cov5_exon <- read.delim("20200129_CpG_cov5_Exon.txt",sep="\t",header = FALSE)
cov5_intron <- read.delim("20200129_CpG_cov5_Intron.txt",sep="\t",header = FALSE)
cov5_gene <- read.delim("20200129_CpG_cov5_gene.txt",sep="\t",header = FALSE)
cov5_intergenic <- read.delim("20200129_CpG_cov5_Intergenic.txt",sep="\t",header = FALSE)
# DML by Trt
DMLTrt_exon <- read.delim("20200129_DML_all_Exon.txt",sep="\t",header = FALSE)
DMLTrt_intron <- read.delim("20200129_DML_all_Intron.txt",sep="\t",header = FALSE)
DMLTrt_gene <- read.delim("20200129_DML_all_gene.txt",sep="\t",header = FALSE)
DMLTrt_intergenic <- read.delim("20200129_DML_all_Intergenic.txt",sep="\t",header = FALSE)
# 0 here
# DML by TP9
DML9_exon <- read.delim("20200129_DML_tp9_Exon.txt",sep="\t",header = FALSE)
DML9_intron <- read.delim("20200129_DML_tp9_Intron.txt",sep="\t",header = FALSE)
DML9_gene <- read.delim("20200129_DML_tp9_gene.txt",sep="\t",header = FALSE)
DML9_intergenic <- read.delim("20200129_DML_tp9_Intergenic.txt",sep="\t",header = FALSE)
# DML by TP80
DML80_exon <- read.delim("20200129_DML_tp80_Exon.txt",sep="\t",header = FALSE)
DML80_intron <- read.delim("20200129_DML_tp80_Intron.txt",sep="\t",header = FALSE)
DML80_gene <- read.delim("20200129_DML_tp80_gene.txt",sep="\t",header = FALSE)
DML80_intergenic <- read.delim("20200129_DML_tp80_Intergenic.txt",sep="\t",header = FALSE)

#### CpG By Feature Summary #####
featureSum <- list()
# All
featureSum$allCG <-c(length(unique(paste0(all_exon$V4,all_exon$V5,all_exon$V6))),
                     length(unique(paste0(all_gene$V4,all_gene$V5,all_gene$V6))),
                     length(unique(paste0(all_intron$V4,all_intron$V2,all_intron$V6))),
                     length(unique(paste0(all_intergenic$V4,all_intergenic$V5,all_intergenic$V6))))
# 5x coverage
featureSum$CG5 <-c(length(unique(paste0(cov5_exon$V4,cov5_exon$V5,cov5_exon$V6))),
                   length(unique(paste0(cov5_gene$V4,cov5_gene$V5,cov5_gene$V6))),
                   length(unique(paste0(cov5_intron$V4,cov5_intron$V2,cov5_intron$V6))),
                   length(unique(paste0(cov5_intergenic$V4,cov5_intergenic$V5,cov5_intergenic$V6))))
# DML Trt
featureSum$DMLTrt <-c(length(unique(paste0(DMLTrt_exon$V4,DMLTrt_exon$V5,DMLTrt_exon$V6))),
                      length(unique(paste0(DMLTrt_gene$V4,DMLTrt_gene$V5,DMLTrt_gene$V6))),
                      length(unique(paste0(DMLTrt_intron$V4,DMLTrt_intron$V2,DMLTrt_intron$V6))),0)
# Note intergenic region is 0 because file is empty so no cases here
# DML Trt x Tp9
featureSum$DML9 <- c(length(unique(paste0(DML9_exon$V4,DML9_exon$V5,DML9_exon$V6))),
                     length(unique(paste0(DML9_gene$V4,DML9_gene$V5,DML9_gene$V6))),
                     length(unique(paste0(DML9_intron$V4,DML9_intron$V2,DML9_intron$V6))),
                     length(unique(paste0(DML9_intergenic$V4,DML9_intergenic$V5,DML9_intergenic$V6))))
# DML Trt x Tp80
featureSum$DML80 <-c(length(unique(paste0(DML80_exon$V4,DML80_exon$V5,DML80_exon$V6))),
                     length(unique(paste0(DML80_gene$V4,DML80_gene$V5,DML80_gene$V6))),
                     length(unique(paste0(DML80_intron$V4,DML80_intron$V2,DML80_intron$V6))),
                     length(unique(paste0(DML80_intergenic$V4,DML80_intergenic$V5,DML80_intergenic$V6))))

Categories <- c("CpGs_all","CpGs_5x","DML_trt","DML_trt_9","DML_trt_80")
Feature <- c("Exon","Gene","Intron","Intergenic")
cat_vec <- NULL
feat_vec <- NULL
count_vec <- NULL
percent_vec <- NULL
for(i in 1:length(Categories)){
  for(j in 1:length(Feature)) {
    cat_vec <- c(cat_vec,Categories[i])
    feat_vec <- c(feat_vec,Feature[j])
    count_vec <- c(count_vec,featureSum[[i]][j])
    percent_vec <- c(percent_vec,round(featureSum[[i]][j]/sum(featureSum[[i]][c(1,3,4)])*100,digits = 2))
  }
}
sumtable <- data.frame(category=cat_vec,feature=feat_vec,count=count_vec,percent=percent_vec)
# Save table
setwd(outputDir)
write.csv(sumtable,"20200129_AllCpGAmongFeaturesSummaryTable.csv")

