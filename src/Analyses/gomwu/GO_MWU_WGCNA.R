
#### Script for running GO-MWU with WGCNA membership values ####

## Here we look at GO enrichment of the 4 modules that were significantly associated with 
 # delta pH (cyan, lavenderblush3,darkred,lightgreen). 

## Description from github repo:
# GO_MWU uses continuous measure of significance (such as fold-change or -log(p-value) ) 
# to identify GO categories that are significantly enriches with either up- or down-regulated genes. 
# The advantage - no need to impose arbitrary significance cutoff. If the measure is binary (0 or 1) 
# the script will perform a typical "GO enrichment" analysis based Fisher's exact test: it will show 
# GO categories over-represented among the genes that have 1 as their measure. On the plot, different 
# fonts are used to indicate significance and color indicates enrichment with either up (red) or down 
# (blue) regulated genes. No colors are shown for binary measure analysis. The tree on the plot is 
# hierarchical clustering of GO categories based on shared genes. Categories with no branch length 
# between them are subsets of each other. The fraction next to GO category name indicates the fracton 
# of "good" genes in it; "good" genes being the ones exceeding the arbitrary absValue cutoff (option 
# in gomwuPlot). For Fisher's based test, specify absValue=0.5. This value does not affect statistics 
# and is used for plotting only.
# Mikhail V. Matz, UT Austin, February 2015; matz@utexas.edu


### Input variables
  # The ModuleMembership .csv files were generated in the '06_CreatingWGCNAObj.R' script in the RNAseq data folder.
  # Path to previous script: 'AE17_Cvirginica_MolecularResponse/src/RNA_seq/06_CreatingWGCNAObj.R'
  # Path to ModuleMembership files (csv) : 'AE17_Cvirginica_MolecularResponse/data/Analysis'
  # Path to goAnnotations  and goDatabase : 'AE17_Cvirginica_MolecularResponse/data/Analysis/gomwu'

## Note to rerun this you will need to move files into the same folder or specify unique paths for each variable.

input = list.files(pattern="RNA_Limma_WGCNA_ModuleMembership")[1:4]
goAnnotations=paste0("GOMWU_LogFoldChangeGE_Goterms.tab")
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision = c("MF","BP","CC") # either MF, or BP, or CC
source("gomwu.functions.R")

# Extracts module names from file names
names <- NULL
for(i in 1:length(input)){
      names <- c(names,substring(input[i],34,nchar(input[i])-4))  
}

# Perform the gomwuStats() function on all modules for each go category
# Couple of notes on arguements:
    
    # We use a `clusterCutHeight=0.25` which means that a group of 
    # categories will be merged if the most dissimilar two of them 
    # share >75% of genes included in the smaller of the two
    
    # We specify the 'module=TRUE' which is required when looking at
    # WGCNA module membership ("kME") outputs. In our case we only 
    # use the arguement 'Module=TRUE' and NOT 'Alternative="g"' because
    # the kME values were transformed in an earlier step to be signless.

for(i in 1:length(input)){
  for(j in goDivision){
        if( i != 3 && j != "BP") {
        gomwuStats(input[i], goDatabase, goAnnotations,j,
                perlPath="perl",
                largest=0.1, # a GO category will not be considered if it contains more than this fraction of the total number of genes 
                smallest=5, # a GO category should contain at least this many genes to be considered
                Module=TRUE, # This argument is used because we are looking at WGCNA outputs
                clusterCutHeight=0.25) # threshold for merging similar (gene-sharing) terms. See README for details.)
        }
  }
}

## Note there were not enough observations in "lavenderblush3" - "BP" for it to run properly. AS a result this
 # was skipped in the loop.

## Loop for creating gomwu plots
# for(i in 1:length(input)){
#   for(j in goDivision){
#         pdf(paste0("GOMWU_WGCNA_Membership_",names[i],"_",j,".pdf"))
#         gomwuPlot(input[i],goAnnotations,j,
#                   absValue=2,  # genes with the measure value exceeding this will be counted as "good genes".
#                   level1=0.1, # FDR threshold for plotting.
#                   #Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
#                   level2=0.05, # FDR cutoff to print in regular (not italic) font.
#                   level3=0.01, # FDR cutoff to print in large bold font.
#                   txtsize=1.2, # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
#                   treeHeight=0.5 # height of the hierarchical clustering tree
#                   # colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
#         )
#         dev.off()
#     }
# }

## Script for taking outputs from gomwu and combining them into a single csv output.

# Read in GO-MWU outputs from gomwuStats()
fo_BP <- list.files(pattern="MWU_BP_RNA")
fo_MF <- list.files(pattern="MWU_MF_RNA")
fo_CC <- list.files(pattern="MWU_CC_RNA")
fo <- c(fo_MF,fo_BP,fo_CC)

# Take module names from the file name
names <- NULL
for(i in 1:length(fo)){
      names <- c(names,substring(fo[i],41,nchar(fo[i])-4))  
}
# Take the GO category from the file name
cate <- NULL
for(i in 1:length(fo)){
      cate <- c(cate,substring(fo[i],5,6))
}
# Loop through files and generate single large data.frame
for(i in 1:length(fo)){
    temp <- read.delim(fo[i],header=TRUE,sep=" ")
    if(i == 1){comb <- data.frame(Module=names[i],Category=cate[i],temp)}
    if(i>1){comb <- rbind(comb,data.frame(Module=names[i],Category=cate[i],temp))}
}
# Only keep those that are significant
comb <- comb[comb$p.adj <= 0.05,]
# Save file
write.csv(comb,"GOMWU_WGCNA_finalSummary.csv")