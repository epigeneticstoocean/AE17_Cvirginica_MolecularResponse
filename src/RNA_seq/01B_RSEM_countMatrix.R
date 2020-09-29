# Script intended for making a count matrix from the single sample "isoform.results" and "gene.results"
# outputs created by RSEM. These represent quantified transcript counts based on '.bam' files 
# created from a 2Pass STAR mode.

# OPERATION: Given the standard 2PASS STAR followed by RSEM transcript quantification
# there should be a new folder with the RSEM ouputs including the '*isoform.results' and
# '*gene.results'. The only thing that should be done is setting the working directore (i.e. setwd())
# to the appropriate folder with these outputs. 

# OUTPUT: the script will automatically generate outputs for both isoform ('RSEM_isoform_CountMatrix')
# and gene ('RSEM_gene_CountMatrix') level quantification. The script will produce three '.csv' of
# the count matrix, repsented as either estimated counts ('estimatedCount'), transcripts per million
# ('TPM'), or fragments per kilobase ('FPKM'). In addition, a list will be created for each quantification
# level (gene or isoform) and saved as an '.RData' file.
# file, which contains the two counts matrix versions, a two vectors of gene or isoform lengths
# (actual and effective), for isoforms a isoform-gene link matrix (2 colums, first with isoform IDs,
# second with the gene ID),  

# Structure of Count Matrix
#	- ROWS = # of genes or isoforms
#	- COLS = # of samples 

# Structure of '.RData'List
# 	- $Count = Matrix of estimated counts
# 	- $TPM = Matrix of TPM
#	- $FPKM = Matrix of FPKM
#	- $length = length of genes (bp)
#	- $eff_length = Effective length of genes
#	- $ID_link = (2x#oftranscripts) matrix with link between transcript ID and gene_ID
#	- $IsoPct = Matrix of Isoform percentage(only for isoform files)

# This will need to change if rerun for another dataset
setwd("/shared_lab/20180226_RNAseq_2017OAExp/RNA/STAR_files/run_20190726_RSEMGnomonIndex/RSEM")

# Get list of all files (gene and isoforms)
fl_genes <- list.files(pattern="*genes.results")
fl_iso <- list.files(pattern="*isoforms.results")

# Create vector of sample names
file_sep <- matrix(unlist(strsplit(fl_genes,split="_")),ncol=2,byrow=TRUE)

# For Genes
temp <- read.delim(fl_genes[1],header=TRUE)
count_mat <- matrix(0,nrow=nrow(temp),ncol=length(fl_genes))
TPM_mat <- matrix(0,nrow=nrow(temp),ncol=length(fl_genes))
FPKM_mat <- matrix(0,nrow=nrow(temp),ncol=length(fl_genes))
for(i in 1:length(fl_genes)){
	temp_file <- read.delim(fl_genes[i],header=TRUE)
	count_mat[,i] <- temp_file$expected_count
	TPM_mat[,i] <- temp_file$TPM
	FPKM_mat[,i] <- temp_file$FPKM
}
count_df <- data.frame(count_mat)
TPM_df <- data.frame(TPM_mat)
FPKM_df <- data.frame(FPKM_mat)

rownames(count_df) <- temp[,1]
colnames(count_df) <- file_sep[,1]
rownames(TPM_df) <- temp[,1]
colnames(TPM_df) <- file_sep[,1]
rownames(FPKM_df) <- temp[,1]
colnames(FPKM_df) <- file_sep[,1]

write.csv(count_df,"RSEM_gene_EstCount.csv")
write.csv(TPM_df,"RSEM_gene_TPM.csv")
write.csv(FPKM_df,"RSEM_gene_FPKM.csv")

gene_list <- list(Count=count_df,TPM=TPM_df,FPKM=FPKM_df,Length=temp$length,
		  Eff_length=temp$effective_length,ID_Link=temp[,1:2])
saveRDS(gene_list,"RSEM_gene_Summary.Rdata")


# For Isoforms
temp <- read.delim(fl_iso[1],header=TRUE)
count_mat <- matrix(0,nrow=nrow(temp),ncol=length(fl_iso))
TPM_mat <- matrix(0,nrow=nrow(temp),ncol=length(fl_iso))
FPKM_mat <- matrix(0,nrow=nrow(temp),ncol=length(fl_iso))
IsoPct_mat <- matrix(0,nrow=nrow(temp),ncol=length(fl_iso))
for(i in 1:length(fl_iso)){
        temp_file <- read.delim(fl_iso[i],header=TRUE)
        count_mat[,i] <- temp_file$expected_count
        TPM_mat[,i] <- temp_file$TPM
        FPKM_mat[,i] <- temp_file$FPKM
	IsoPct_mat[,i] <- temp_file$IsoPct
}
count_df <- data.frame(count_mat)
TPM_df <- data.frame(TPM_mat)
FPKM_df <- data.frame(FPKM_mat)
IsoPct_df <- data.frame(IsoPct_mat)

rownames(count_df) <- temp[,1]
colnames(count_df) <- file_sep[,1]
rownames(TPM_df) <- temp[,1]
colnames(TPM_df) <- file_sep[,1]
rownames(FPKM_df) <- temp[,1]
colnames(FPKM_df) <- file_sep[,1]
rownames(IsoPct_df) <- temp[,1]
colnames(IsoPct_df) <- file_sep[,1]

write.csv(count_df,"RSEM_isoform_EstCount.csv")
write.csv(TPM_df,"RSEM_isoform_TPM.csv")
write.csv(FPKM_df,"RSEM_isoform_FPKM.csv")

isoform_list <- list(Count=count_df,TPM=TPM_df,FPKM=FPKM_df,IsoPct=IsoPct_df,
	Length=temp$length,Eff_length=temp$effective_length,ID_Link=temp[,1:2])
saveRDS(isoform_list,"RSEM_isoform_Summary.Rdata")

