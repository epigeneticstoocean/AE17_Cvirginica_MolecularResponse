library(data.table)
library(dplyr)

sprintf("Reading in data and initializing methylation count matrix...")
# Set working directory to base folder where your raw files and bismark outputs reside
setwd("/shared_lab/20180226_RNAseq_2017OAExp/DNAm/trimmedSamples/singleSample_trimScript_20190719")

# Specific folder where cov matrix are.
bis_outputs <- "CytoReports"
file_outputs <- "countMatrix"

files <- list.files(bis_outputs)
count_files <- files[grep("CytoSummary.CpG_report.txt",files)]
samples <- substr(count_files,1,5)

init <- fread(paste0(bis_outputs,"/",samples[1],"_CytoSummary.CpG_report.txt"),header=FALSE)

sum.table <- init[,c(1,2,3,6,7)]

init_C <- init[,4]
colnames(init_C) <- samples[1]
init_T <- init[,5]
colnames(init_T) <- samples[1]

sprintf("Starting counting....")
n<- 23
pb <- txtProgressBar(min = 0, max = n, style = 3)

for(i in 2:length(samples)){
  Sys.sleep(0.001)
  setTxtProgressBar(pb, i-1)

  temp <- fread(paste0(bis_outputs,"/",samples[i],"_CytoSummary.CpG_report.txt"),header=FALSE)
  temp_C <- temp[,4]
  temp_T <- temp[,5]

  colnames(temp_C) <- samples[i]
  colnames(temp_T) <- samples[i]

  init_C <- cbind(init_C,temp_C)
  init_T <- cbind(init_T,temp_T)
}
close(pb)

sprintf("Completed merging and reformatting matrices...")

#Calculating Methylation
final_S <- init_C + init_T
final_B<- init_C/final_S

sprintf("Saving files...")
fwrite(sum.table,paste0(file_outputs,"/All_CytoSum_summaryTable.csv"))
saveRDS(sum.table,paste0(file_outputs,"/All_CytoSum_summaryTable.RData"))


fwrite(init_C,paste0(file_outputs,"/All_CytoSum_methylCountMatrix.csv"))
fwrite(init_T,paste0(file_outputs,"/All_CytoSum_unmethylCountMatrix.csv"))
fwrite(final_S,paste0(file_outputs,"/All_CytoSum_TotalCountMatrix.csv"))
fwrite(final_B,paste0(file_outputs,"/All_CytoSum_BetaMatrix.csv"))

saveRDS(init_C,paste0(file_outputs,"/All_CytoSum_methylCountMatrix.RData"))
saveRDS(init_T,paste0(file_outputs,"/All_CytoSum_unmethylCountMatrix.RData"))
saveRDS(final_S,paste0(file_outputs,"/All_CytoSum_TotalCountMatrix.RData"))
saveRDS(final_B,paste0(file_outputs,"/All_CytoSum_BetaMatrix.RData"))
