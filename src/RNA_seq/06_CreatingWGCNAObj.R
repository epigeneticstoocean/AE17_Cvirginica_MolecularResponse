
## Analysis for performing weighted gene co-expression network analysis
 # The objective of this script is to use the WGCNA approach to cluster genes that
 # co-express in order to find associations between environmental and phenotypic 
 # variables and expression data when gene expression data is subtle.
 # Outputs : ...


##### Section Run on Cluster ######

## Library
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
library(limma)
library(cowplot)


#### Data ####
# Set wd()
# Display the current working directory
getwd()
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = "/home/downeyam/Github/AE17_Cvirginica_MolecularResponse"
#workingDir = "."
setwd(workingDir)

# Meta data
traitData <- readRDS("data/meta/AE17_RNAmetaData.RData")
traitData_sans17005 <-  traitData[traitData$ID != "17005",]
dim(traitData_sans17005)
names(traitData_sans17005)
#Read in the gene expression data
dataExp <- readRDS("results/RNA/RNA_gene_postVoomAndNormalization_DGEListObj.RData")
gc <- dataExp$E
gc_sans17005 <- gc[,traitData$ID != "17005"]
# Removing 17005 outlier
datExpr0 <-  as.data.frame(t(gc_sans17005))
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK
#[1] TRUE

#### Trying a couple of different clustering approaches ####
sampleTree = hclust(dist(datExpr0), method = "average")
sampleTree1 = hclust(dist(datExpr0), method = "centroid")
sampleTree2 = hclust(dist(datExpr0), method = "ward.D2")
sampleTree3 = hclust(dist(datExpr0), method = "median")
# Storing them as a list
trees <- list(sampleTree,sampleTree1,
              sampleTree2,sampleTree3)
# Vector of the different clustering approaches
approaches <- c("Average","Centroid","Ward.D2","Median")

# Create plots for each clustering approach
for(i in 1:length(trees)){
        pdf(file = paste0("WGCNA_cluster_",approaches[i],"method.pdf"), width = 12, height = 9)
        plot(trees[[i]],main=approaches[i])
        dev.off()
}
# Looks like using the WARD.D2 clustering algorithm performs the best so I'll proceed with that
# heirarchical clustering strategy for all downstream steps.

#Determine cluster under the line
# clust = cutreeStatic(sampleTree, cutHeight = 200, minSize = 10)
# table(clust)
datExpr <- datExpr0
# Creating trait object from metadata file
traitD <- traitData_sans17005[,c(6:9,20:30)]
traitD[,11:15] <- sapply(traitD[,11:15],as.character)
traitD[,11:15] <- sapply(traitD[,11:15],as.numeric)
rownames(traitD) = traitData_sans17005$ID

# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(traitD, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
pdf(file = "Dendrogram.pdf", width = 12, height = 9)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(traitD), 
                    main = "Sample dendrogram and trait heatmap")
dev.off()
save(datExpr,traitD, file = "results/RNA/RNA_Limma_Expression_Data_forWGCNA.RData")

#### Step 2 : Clustering ####
# Make sure to set working directory were where the RData file is located
library(WGCNA)
enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "results/RNA/RNA_Limma_Expression_Data_forWGCNA.RData");
# The variable lnames contains the names of loaded variables.
lnames

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
pdf("softThresholdPowerIndice2.pdf")
#sizeGrWindow(9, 5)
#par(mfrow = c(1,2));
#cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(-sign(sft$fitIndices[,3])*sft$fitIndices[,2]~sft$fitIndices[,1],type="l",
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",
     main = paste("Scale independence"))

# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
dev.off()

png("meanConnectivityAndPower.png")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,col="red")
dev.off()

# This soft power was choose because it was the
# threshold were 90% of data was explained
softPower = 5;
adjacency = adjacency(datExpr, power = softPower);
# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "ward.D2")
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
pdf("GeneCluster_TOM_ward2.pdf")
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
dev.off()

geneTree2 = hclust(as.dist(dissTOM), method = "average")
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
pdf("GeneCluster_TOM_average.pdf")
plot(geneTree2, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
dev.off()

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
pdf("Dendrogram.pdf")
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "ward.D2");
# Plot the result
sizeGrWindow(7, 6)
pdf("METree_Clustering_ward.pdf")
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "results/RNA/RNA_Limma_networkConstruction_WGCNA.RData")

#### Extract Module Membership values for GO-MWU ####

lnames = load(file = "results/RNA/RNA_Limma_Expression_Data_forWGCNA.RData")
lnames = load(file = "results/RNA/RNA_Limma_networkConstruction_WGCNA.RData")
MEs_all = moduleEigengenes(datExpr, moduleColors)
nSamples = nrow(datExpr)
datME <- moduleEigengenes(datExpr,moduleColors)$eigengenes
geneList <- colnames(datExpr)
MEs = orderMEs(MEs0)
# Calculate correlation among target variables and modules
moduleTraitCor = cor(MEs, traitD, use = "p")
modTrait_corr <- data.frame(moduleTraitCor)
# Calculate P value from correlation
modTrait_P <- data.frame(corPvalueStudent(moduleTraitCor, nSamples))
# List of module gene summary information
modList <- list(datME,moduleColors,modTrait_P,modTrait_corr)
topDiffpHNames <- rownames(modTrait_P)[order(modTrait_P$diff_pH)][2:5]

MEs_target <- MEs[,colnames(MEs) %in% topDiffpHNames]

out <- signedKME(datExpr,
                 MEs_target,
                 outputColumnName = "",
                 corFnc = "cor",
                 corOptions = "use ='p'")

moduleMemberList <- list()
for( i in 1:length(topDiffpHNames)){
        gl <- geneList[paste0("ME",modList[[2]]) == topDiffpHNames[i]]
        temp <- data.frame(gene_id = rownames(out),
                           Mod_Membership = out[,colnames(out) %in% substring(topDiffpHNames[i],3)])
        temp$Mod_Membership[!temp$gene_id %in% gl] <- 0
        temp$Mod_Membership[temp$gene_id %in% gl] <- abs(temp$Mod_Membership[temp$gene_id %in% gl])
        moduleMemberList[[i]] <- temp
        write.csv(temp,paste0("results/RNA/RNA_Limma_WGCNA_ModuleMembership_",
                         substring(topDiffpHNames[i],3),".csv"),row.names = FALSE)
}
names(moduleMemberList) <- substring(topDiffpHNames,3)
        
saveRDS(moduleMemberList,"results/RNA/RNA_Limma_WGCNA_ModuleMembership.RData")
