library(parallel)

setwd("/shared_lab/20180226_RNAseq_2017OAExp/DNAm")
#source("scripts/R/DNAmRefCode.R")

sprintf("Reading in Data....")
# Gene information
gene <- readRDS("reference/gene_GeneLoc.RData")
# Create new gene data.frame that has additional index column based on row names
y <- data.frame(index=1:nrow(gene),gene)

# CpG meta data (based on cytosummary flag in bismark)
meta <- readRDS("processed_samples/04_countMatrix/All_CytoSum_summaryTable.RData")
# Create a unique index integer for each cytosine in CytoSummary File
meta$index <- c(1:nrow(meta))

# Set cores for running in parallel
ncores <- unique(y$chr)

# Loop through each gene and find CpGs that
# fall between the gene range.
# Data passed:
# CpG info  : chr, position, and strand 
# Gene info : start , end

findMatch <- function(index,chr,m,f,gID,gid,strand,start,end){
    if(any(which(as.numeric(tmp_cpg$V2) >= start))){
      temp <- tmp_cpg[which(as.numeric(tmp_cpg$V2) >= start),]
      if(any(which(as.character(temp$V3) == strand))){
        temp <- temp[which(as.character(temp$V3) == strand),]
        if(any(which(as.numeric(temp$V2) <= end))){
          temp <- temp[which(as.numeric(temp$V2) <= end),]
	  #out_temp<-data.frame(index,chr,m,f,gID,gid,
          #                     strand,start,end,
	  #		       order=c(1:nrow(temp)),temp)
	  out_temp<-cbind(as.character(index),chr,m,f,gID,gid,
                               strand,start,end,
                               order=c(1:nrow(temp)),temp)
	  return(out_temp)
        }
      }
    }
}

mapply_findMatch <- function(id,meta){
  # Subset of genes for a particular chromosome
  tmp_gene <- y[id[1]:id[2],]
  #Only process cpgs within the chromosome for faster processing time
  tmp_cpg <- meta[meta$V1 == unique(tmp_gene$chr),]
  final<-mapply(findMatch,tmp_gene$index,tmp_gene$chr,tmp_gene$method,
		   tmp_gene$feature,tmp_gene$gene_ID,tmp_gene$gene_id,
		   tmp_gene$strand,tmp_gene$start,tmp_gene$end)
  return(final)
}

idm<-matrix(ncol=2,nrow=length(unique(y$chr)))
for(i in 1:length(unique(y$chr))){
	temp<-y[y$chr == unique(y$chr)[i],]
	idm[i,] <- cbind(min(temp$index),max(temp$index))
}

sprintf("Starting matching.....")
sprintf(paste0("Working on..",ncores," cores"))

result <-mclapply(nrow(idm):1,
                    function(x) mapply_findMatch(idm[x,]),
                    mc.cores=ncores)
sprintf("Preliminary Index List Saved...")
saveRDS(result,"CytoSummary_Gene_Index_RAWList.RData")

#sprintf("Removing nulls...")
#result <- mclapply(result,
#function(x){x[x!="NA"]},
#function(x){Filter(Negate(is.null),x)},
#mc.cores=ncores)

sprintf("Unlisting and creating matrix...")
index <- matrix(unlist(result),ncol=3,byrow = TRUE)
saveRDS(index,"CytoSummary_Gene_IndexwithNA.RData")

sprintf("Condensing index list (removing NAs)...")
index <- index[index[,1] != "NA",]
saveRDS(index,"CytoSummary_Gene_Index.RData")

sprintf("Subsetting  meta data and gene list with indexes...")
full <- cbind(meta[index[,1],],index[,2],gene[index[,3],])
col<-c("chr_1","pos","cg_strand","motif","tri",
"cg_index","match_num","gene_index","chr_2","method",
"feature","start","end","score","gene_strand","phase","gene_ID",
"gene_id")
colnames(full)<-col

sprintf("Saving final data set....")
saveRDS(full,"CytoSummary_Gene_Match.RData")

