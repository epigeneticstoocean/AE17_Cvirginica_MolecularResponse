#### This script is used summarize DNA Methylation response to OA and time ####
 # Figures from this code : 3

## Libraries
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(data.table)
library(vegan)
library(adegenet)
library(matrixStats)
# color palette
pal <- brewer.pal(n = 12, name = 'Paired')
col_perm <- c(pal[1],pal[5],pal[2],pal[6])

#### Data ####
##Path
inputDir <- "~/Github/AE17_Cvirginica_MolecularResponse/"
setwd(inputDir)

### Sample Meta Data ###  
meta <- readRDS("data/meta/AE17_RNAmetaData.RData")
meta <- meta[meta$ID != "17099",] # need to remove individual 17099

### Summary table of Cpgs by feature (Panel 1)###
 # Contains counts of CpGs by the different primary subsets (all,those with coverage of >= 5,
 # those from the different DML comparisons) among exons,introns, and intergenic regions (currently didn't 
 # do anything more than those groups)
sumTable <- read.csv("results/DNAm/DNAm_20200129_AllCpGAmongFeaturesSummaryTable.csv")
sumTable_rev <- sumTable[sumTable$feature != "Gene",] # Remove gene in this case because it the same as introns + exons
# Restructuring the table for plotting
sumTable_rev$feature <- factor(sumTable_rev$feature,levels=c("Exon","Intron","Intergenic" ))
sumTable_rev$category <- factor(sumTable_rev$category,levels=c("CpGs_all","CpGs_5x","DML_trt","DML_trt_9","DML_trt_80"))
levels(sumTable_rev$category) <- c("CpGs (all)","CpGs (5x)",
                                   "DML","DML (D9)","DML (D80)")

### All CpGs by feature ###
# Example: meth$mC$exon would retrieve a table for
#          the methylated cytosine counts (mC) within exons 
#meth <- readRDS("results/DNAm/20200202_AllCountsList_cov5_byFeature.RData")
#m <- meth$beta
#saveRDS(m,"results/DNAm/20200202_AllBeta_cov5_byFeature.RData")
# # total Count data for genes
# meth_total_gene <- meth$tC$gene[,2:24] # remove first column which has coordinate information
# sumCountMethylation_gene <- colSums(meth_total_gene,na.rm = TRUE)
# # beta values for genes
#meth_beta <- readRDS("results/DNAm/20200202_AllBeta_cov5_byFeature.RData")
#meth_beta_gene <- meth_beta$gene[,2:24]
meth_beta_gene <- readRDS("results/DNAm/20200202_geneBeta_cov5_byFeature.RData")
### Calculate median methylation globally for each individual for each feature
med_values <- NULL
name_values <- NULL
for(i in 1:length(meth_beta)){
  temp <-  as.matrix(meth_beta[[i]][,2:24])
  nm <- rep(names(meth_beta)[i],times=23)
  temp_medians <- colMedians(temp,na.rm = TRUE)
  
  med_values <- c(med_values,temp_medians)
  name_values <- c(name_values,nm)
}
IDs <- rep(meta$ID,length(meth_beta))
Treatment <- rep(meta$Treatment,length(meth_beta))
Time <- rep(meta$Time,length(meth_beta))
median_df <- data.frame(IDs,Treatment,Time,feature=name_values,median=med_values)
# Remove gene since it is redundant (exon + intron)
median_dfnoGene <- median_df[median_df$feature != "gene",]
median_dfnoGene$feature <- as.factor(as.character(median_dfnoGene$feature))
# Reorder levels of feature
median_dfnoGene$feature <- factor(median_dfnoGene$feature,levels=c("exon","Intron","Intergenic" ))
levels(median_dfnoGene$feature) <- c("Exon","Intron","Intergenic")

methByGene <-  readRDS("data/MBDBS_seq/20200130_CpGbyGeneSummary/gene_CpGcoverageSummary.RData")

#### Figure Three ####

#### Plot 1 - CpG by feature summary across CpG subsets ####

P1 <- ggplot(sumTable_rev,aes(x=category,y=percent,fill=feature)) + geom_bar(stat = "identity") + 
  theme_cowplot() + scale_fill_manual(values = viridisLite::viridis(3)) +
  labs(y="CpGs by Feature (%)",x="",title="\n\n",fill="") + 
  ylim(0,100) +
  scale_x_discrete(labels = c(bquote(CpGs[all]),
                              bquote(CpGs[5]),
                              bquote(DML[Trt]),
                              bquote(DML[Trt_D9]),
                              bquote(DML[Trt_D80]))) +
  labs(x = "",y="CpGs by Feature (%)") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size=14, angle=45, hjust=1, vjust=1))
P1
## NOTE Summary of the total number of CpGs in each column (subset) was done 
 # manually in inkscape.
sumTableSub <- sumTable[sumTable$feature == "Gene" | sumTable$feature == "Intergenic",]
sumTableSub$count[seq(1,10,by=2)]+ sumTableSub$count[seq(2,10,by=2)]
# Number of CpGs for each column of panel 1

#### Plot 2 - Global Median Methylation among TreatmentsXTimes ####
P2 <- ggplot(median_dfnoGene,
             aes(x=interaction(Treatment,Time),
                 y=median*100)) + 
  geom_boxplot(aes(fill = interaction(Treatment,Time))) + 
  theme_cowplot() + 
  facet_grid(.~feature) + 
  ylim(75,95) +
  scale_fill_manual(values = col_perm) + 
  labs(x = "",y="Median Methylation (%)") + 
  theme(legend.position = "none",
        axis.title.x=element_blank(), # Comment these lines if you want x labels
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
P2

# Full model with feature
model_out <- aov(median~Treatment*Time*feature,data=median_dfnoGene)
(out <- summary(model_out))

#### Plot 3 -  PCA plots separated by feature ####
# 1 = Exon, 2 = gene, 3 = Intron, 4 = intergenic region
pca_plots <- list()
plot_titles <- c("Exons","Gene Body","Introns","Intergenic")
for(i in 1:length(meth_beta)){
  temp <- as.matrix(meth_beta[[i]][,2:ncol(meth_beta[[i]])])
  temp[is.na(temp)] <- 0
  class(temp) <- "numeric"
  pca <- prcomp(t(temp))
  eigs <- pca$sdev^2
  temp_df <- data.frame(x=pca$x[,1],y=pca$x[,2],condition=meta$SFV)
  p_temp <- ggplot(temp_df,aes(x=x,y=y,colour=condition,shape=condition)) + 
    geom_point(size=4) + 
    theme_cowplot() + 
    scale_colour_manual(values = col_perm[c(1,3,2,4)],
                        labels=c('Control\n     D9','Control\n     D80','High OA\n     D9', 'High OA\n     D80')) + 
    scale_shape_manual(values = c(15,15,17,17),
                       labels=c('Control\n     D9','Control\n     D80',
                                'High OA\n     D9', 'High OA\n     D80')) +
    labs(x=paste0("PC1 (",round(eigs[1] / sum(eigs)*100,1),"%)"),
         y=paste0("PC2 (",round(eigs[2] / sum(eigs)*100,1),"%)"),
         title = plot_titles[i],
         colour="",
         shape="") +
    theme(plot.title = element_text(hjust = 0.5))
  pca_plots[[i]] <- p_temp + theme(legend.position="none") 
}
# Extract the legend from one of the PCA plots
legend <- get_legend(
 # create some space to the left of the legend
 p_temp + theme(legend.box.margin = margin(0, 0, 0, 12),
                legend.title = element_blank())
)

## Significant using adonis function from vegan for gene bodies
gb_dnam <- temp <- as.matrix(meth_beta[[2]][,2:ncol(meth_beta[[2]])])
gb_dnam[is.na( gb_dnam)] <- 0
class( gb_dnam) <- "numeric"
(out_dnam <- adonis(t(gb_dnam)~Time:Treatment+Time+Treatment+Pop+Lane,data=meta,
                  permutations = 5000,method = "manhattan" ))
# Subtle effect of treatment but not effect of time or the interaction

#### Plot 4 - DAPC separated by feature ####
dapc_plot <- function(x,model_dnam,names=NULL) {
  #x <- meth_beta$gene
  x <- as.matrix(x[,2:ncol(x)])
  x[is.na(x)] <- 0
  class(x) <- "numeric"
  early_time_dnam <- x[,model_dnam$Day == 9]
  early_time_meta_dnam <- model_dnam[model_dnam$Day == 9,]
  dapc_treatment <- dapc(t(early_time_dnam),early_time_meta_dnam$treatment,n.pca=7,n.da=1)
  early_time_meta_dnam$coord<- unlist(dapc_treatment$ind.coord[,1])
  ## Mapping Day 80
  late_time_dnam <- x[,model_dnam$Day == 80]
  late_time_meta_dnam <- model_dnam[model_dnam$Day == 80,]
  predict_values <- predict.dapc(dapc_treatment,t(late_time_dnam))
  late_time_meta_dnam$coord <-unlist(predict_values$ind.scores[,1])
  whole_meta_dnam<- rbind(early_time_meta_dnam,late_time_meta_dnam)
  whole_meta_dnam$coord_transform <- log(whole_meta_dnam$coord+abs(min(whole_meta_dnam$coord))+1)
  
  output <- ggplot(whole_meta_dnam,aes(x=coord,fill=SFV)) + 
            geom_density(adjust=2) + xlim(-8,8) +
            scale_fill_manual(values = col_perm[c(1,3,2,4)],
                                labels=c('Control\n     D9','Control\n     D80',
                                         'High OA\n     D9', 'High OA\n     D80')) +
            scale_y_continuous(breaks = c(0.1,1,2,4,6,10,20),limits = c(0,20),trans = "sqrt") + 
            theme_cowplot() +
            labs(x="Coordinate",
                 y="Density (sqrt)",
                 title = names[i]) +
            theme(plot.title = element_text(hjust = 0.5),
                  legend.title = element_blank()) 
  return(output)
}

dapc_plots <- list()
plot_titles <- c("Exons","Gene Body","Introns","Intergenic")
for(i in 1:length(meth_beta)){
  dapc_plots[[i]] <- dapc_plot(meth_beta[[i]],meta,plot_titles)
}
legend_dapc <- get_legend(
  # create some space to the left of the legend
  dapc_plots[[1]] + theme(legend.box.margin = margin(0, 0, 0, 12),
                 legend.title = element_blank())
)


#### Constructing Final Figure 3 ####
## Option Current ##
# First column of figure 3 
# plot 1 - cpg summary by feature for different subsets
# plot 2 - median methylation by treatment level
first <- plot_grid(P1,P2,labels=c("A","D"),ncol=1)
# Second column of figure 3
second <- plot_grid(pca_plots[[2]],
                    dapc_plots[[2]]+ theme(legend.position="none"),
                    labels=c("B","C"),ncol=1)
# Third column - legend
third <- plot_grid(legend,legend_dapc,ncol=1)
second_third <- plot_grid(second,third,rel_widths = c(3, 1.5))
## Final figure
# First combine first two columns
right_column <- plot_grid(second_third,NULL,labels = c("","E"),ncol=1)
# Then add legend (adjusting for correct widths)
final_1 <- plot_grid(first,right_column,ncol=2)
final_1
#ggsave(plot=final_1,filename = "results/manuscript/figures/Figure3/Fig3_final.png")
#ggsave(plot=final_1,filename = "results/manuscript/figures/Figure3/Fig3_final.pdf")

## Final figure with DML venn diagram created in inkscape

#### Supplemental Figure ####

### Genome wide DNA methlyation patterns with PCA and DAPC
pca_col <- plot_grid(plotlist = pca_plots[c(1,3,2,4)],ncol=1,labels=c("A","C","E","G"))
dapc_col <- plot_grid(dapc_plots[[1]]+ theme(legend.position="none"),
                      dapc_plots[[3]]+ theme(legend.position="none"),
                      dapc_plots[[2]]+ theme(legend.position="none"),
                      dapc_plots[[4]]+ theme(legend.position="none"),
                      ncol=1,labels=c("B","D","F","H"))
final_gw_supp <- plot_grid(pca_col,legend,dapc_col,legend_dapc,rel_widths = c(3,1,3,1),ncol=4)
final_gw_supp
#ggsave(final_gw_supp,filename="results/manuscript/Supp/figures/supp_GW_DNAm.png")
#ggsave(final_gw_supp,filename="results/manuscript/Supp/figures/supp_GW_DNAm.pdf")
