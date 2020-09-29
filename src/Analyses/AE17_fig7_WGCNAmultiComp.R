#### WGCNA and Differential pH association ####

## Figures from code: 7

## This script takes the two RData objects generated from the WGCNA script to visualize
 # and analyze the gene co-expression modules with phenotype, DNA methylation, and enviro-
 # nment. 

## Libraries
library(WGCNA)
library(cowplot)
library(ggplot2)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(data.table)
library(GO.db)
library(dplyr)

## colors 
pal <- brewer.pal(n = 12, name = 'Paired')
col_perm <- c(pal[1:2],pal[5:6],pal[12])

#### Data ####
setwd("/home/downeyam/Github/AE17_Cvirginica_MolecularResponse/")

### Meta data ###
model_original<-read.csv("data/Phenotype/CompletePhenotype_final2020.csv",stringsAsFactors = FALSE)
sampleList<-readRDS("data/meta/metadata_20190811.RData")
model <- sampleList[sampleList$ID != "17005",]
m_final <- model[model$ID != "17099",]
mo <- model_original[model_original$ID %in% m_final$ID,]
m_final$EPF_pH <- mo$EPF_pH_Total
m_final$diff_pH <- mo$EPF_pH_Total-mo$pH_Total_2W

ref <- readRDS("data/references/CDS_wGoTerms_GeneLOC.RData")
ref <- ref[!duplicated(ref$gene_id),]

### Module Gene Expression Data ###
# Files from the WGCNA R script (run on separately on a computing cluster)
# Load the expression and trait data saved in the first part
lnames = load(file = "results/RNA/Limma_Expression_Data_forWGCNA.RData");
#The variable lnames contains the names of loaded variables.
# Load network data saved in the second part.
lnames = load(file = "results/RNA/RNA_Limma_networkConstruction_WGCNA.RData");
# Number of genes
nGenes = ncol(datExpr)
# Number of samples
nSamples = nrow(datExpr)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
geneList <- colnames(datExpr)
MEs = orderMEs(MEs0)
## Creates new object with just the eigenGene Expression values for all modules
datME <- moduleEigengenes(datExpr,moduleColors)$eigengenes
# Calculate correlation among target variables and modules
moduleTraitCor = cor(MEs, traitD, use = "p")
modTrait_corr <- data.frame(moduleTraitCor)
# Calculate P value from correlation
modTrait_P <- data.frame(corPvalueStudent(moduleTraitCor, nSamples))
# List of module gene summary information
modList <- list(datME,moduleColors,modTrait_P,modTrait_corr)

## Selecting Target Module Genes
## Decided to highlight the two modules most strongly associated with diff. pH and
# also the one module that was most strongly associated with Treatment.
# Looking at the top association with diff. pH I found that the pattern was driven 
# by a single individual, so I decided to use the second and third ranked modules instead.
# Selecting 2 and 3 top candidates for diff. pH
modTrait_P[order(modTrait_P$diff_pH),]
topDiffpHNames <- rownames(modTrait_P)[order(modTrait_P$diff_pH)][c(2:4,8)]
# Selecting best candidate for Treatment
#topTreatmentNames <- rownames(modTrait_P)[order(modTrait_P$Treatment)][1]
# Combining targets into single vector
topModuleNames <- c(topDiffpHNames)#,topTreatmentNames)

### DNA MEthylation Data (CpGs) ###
  # NOTE: this file needs to be uncompressed prior to running this line
gene_LOC <- fread("data/references/CpG_5x_methylKit_CDSAnnotation.tab",sep = "\t")
cpg <- readRDS("data/MBDBS_seq/methylKitObj_all_cov5Filtered_united.RData")
cpg_t <- getData(cpg)[,c(seq(5,73,3))]
cpg_methyl <- getData(cpg)[,c(seq(6,73,3))]
cpg_label <- getData(cpg)[,c(1,2)]
# Beta calculated as the number of methylated counts / total number of counts
cpg_beta <- data.frame(cpg_methyl/cpg_t)
# Removing 17099 which had poor CpG coverage and 17005 which was an outlier for
# gene expression
colnames(cpg_beta) <- sampleList$ID[sampleList$ID != "17099"]
cpg_beta$cpg_pos <- as.character(paste0(cpg_label$chr,"_",cpg_label$start))
cpg_beta <- cpg_beta[,sampleList$ID != "17005"]
gene_LOC_simple <- gene_LOC[,c(4,8)]
join_cpg <- left_join(gene_LOC_simple,cpg_beta,by = "cpg_pos")

join_cpg %>% group_by(gene_id) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE) -> cpg_geneSummary
join_cpg %>% group_by(gene_id) %>%
  tally() -> cpg_geneSummary_count
# Create single list with summary beta for each gene and a count for each CpG in each gene
cpgs <- list(cpg_geneSummary,cpg_geneSummary_count)

#### Basic Heatmap (not in figure) ####
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(modTrait_P, 1), ")", sep = "");
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traitD),
               yLabels = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               #textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

#### Heat map with modules order (in figure) ####
# Lets order the heat map by diff_pH (our primary phenotype)
moduleTraitPvalue_orderdiffPH <- modTrait_P[order(modTrait_P[,9]),]
moduleTraitCor_orderdiffPH <- moduleTraitCor[rev(order(moduleTraitCor[,9])),]

# Select only the columns we really care about (in this case the epf measures 
# and the environment conditions)
mod <- moduleTraitCor_orderdiffPH[,c(9,11,12)]
#mod_complete <- moduleTraitCor_orderdiffPH[,c(11,12,8,9,10,1,2,3,4)]
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = mod,
               xLabels = dimnames(mod)[[2]],
               yLabels = dimnames(mod)[[1]], 
               cex.lab = 1.4,
               #ySymbols = names(MEs[rev(order(moduleTraitCor[,2]))]),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               #textMatrix = textMatrix,
               #setStdMargins = FALSE,
               #cex.text = 0.5,
               zlim = c(-1,1),
               main = paste(""))

#### Summarizing modules ####
#### Summary Function ####
moduleSummary <- function(meta,targetMods,mod,cpgs){
  ## Create a list of data.frames that contain module id, sample id, treatxtime information,
  # diff. pH data, and the eigen expression value.
  
  ## Performs a complete summary on the selected modules with the mod_list 
  mod_list <- list()
  # Figures
  mod_P_A <- list() # list of DNA Methylation vs Env. plots
  mod_P_B <- list() # list of EigenGene Expression vs. DNA Mehtylation
  mod_P_C <- list() # list of Diff. pH vs.EigenGene Expression
  # Annotations
  mod_list_full <- list() # list of genes within modules and annotations
  # Statistical Analysis
  mod_list_aov_env_cpg_out <- list() # aov environment vs cpg anova
  mod_list_tukey_env_cpg_out <- list() # tukey environment vs. cpg tukey
  mod_list_lm_cpg_expr_out <- list() # lm DNA methylation vs. Expression
  mod_list_lm_expr_diffpH_out <-  list() # lm Expression vs. DNA methylation
  mod_list_stats_outputs <- list()
  #p_list <- list()
  #corr_list <- list()
  
  for(i in 1:length(targetMods)) {
    #i=3
    ### Module list Summary
    mod_list[[i]] <- data.frame(mod=targetMods[i],
                                samples= meta$ID,
                                SFV=meta$SFV,
                                diff_pH=meta$diff_pH,
                                eigen=mod[[1]][model$ID != "17099",targetMods[i]])
    mod_list[[i]]$SFV <- factor(mod_list[[i]]$SFV,levels=c("09.400","09.2800","80.400","80.2800"))
    
    id <- geneList[paste0("ME",mod[[2]]) == targetMods[i]]
    # Gene IDs
    geneTemp <- data.frame(gene_id=id)
    # Adding DNA methylation counts
    cpg_count_mod <- left_join(geneTemp,cpgs[[2]],by="gene_id")
    geneTemp_wCpG <- data.frame(geneTemp,cpg_count=cpg_count_mod$n)
    
    ### Creating full gene annotation for each model
    mod_list_full[[i]] <- left_join(geneTemp_wCpG,ref,by="gene_id")
    # DNA Methylation Summary
    dna_m_temp <- left_join(geneTemp,cpgs[[1]],by="gene_id")
    # Proportion of CpGs with coverage per module
    mod_list[[i]]$cpg_prop <- sum(is.na(dna_m_temp[,2]))/length(dna_m_temp[,2])
    # Mean methylation for each individual for each module
    unlist(dna_m_temp %>% summarise_if(is.numeric, mean, na.rm = TRUE)) -> mod_list[[i]]$cpg_mean
    
    ### Statistical Analysis
    # Env. vs. DNA Methylation (ANOVA and Tukey)
    temp <- aov(c(mod_list[[i]]$cpg_mean*100)~meta$Treatment*meta$Time)
    temp_out <- summary(temp)
    mod_list_aov_env_cpg_out[[i]] <- temp_out
    
    if(temp_out[[1]]$`Pr(>F)`[3] <= 0.05){
      mod_list_tukey_env_cpg_out[[i]] <- TukeyHSD(temp)
    } else{
      mod_list_tukey_env_cpg_out[[i]] <- "No Interaction"
    }
    # DNA Methylation vs. Expression (linear model)
    mod_list_lm_cpg_expr_out[[i]] <- summary(lm(mod_list[[i]]$eigen~c(mod_list[[i]]$cpg_mean*100)))
    # Expression vs. diff pH (linear model)
    mod_list_lm_expr_diffpH_out[[i]] <- summary(lm(mod_list[[i]]$diff_pH~mod_list[[i]]$eigen))
    
    ### Summarizing statistical outputs
    labels <- c("CpG_Treatment_F","CpG_Time_F","CpG_Interaction_F",
    "CpG_Treatment_P","CpG_Time_P","CpG_Interaction_P",
    "Exp_CpG_R2","Exp_CpG_slope","Exp_CpG_F","Exp_CpG_P",
    "DiffpH_Exp_R2","DiffpH_Exp_slope","DiffpH_Exp_F","DiffpH_Exp_P")
    
    ### Summary of mod Gene Expression association p values and correlations
    #p_list[[i]] <- mod[[3]][match(targetMods[i],rownames(mod[[3]])),c(11,12,9)]
    #corr_list[[i]] <- mod[[3]][match(targetMods[i],rownames(mod[[4]])),c(11,12,9)]
    
    anova_outputs <- c(mod_list_aov_env_cpg_out[[i]][[1]]$`F value`[1:3],
                       mod_list_aov_env_cpg_out[[i]][[1]]$`Pr(>F)`[1:3])
    lm_cpg_outputs <- c(mod_list_lm_cpg_expr_out[[i]]$adj.r.squared,
                        unlist(mod_list_lm_cpg_expr_out[[i]]$coefficients[2,c(1,3,4)]))
    lm_Exp_outputs <- c(mod_list_lm_expr_diffpH_out[[i]]$adj.r.squared,
                             unlist(mod_list_lm_expr_diffpH_out[[i]]$coefficients[2,c(1,3,4)]))
    mod_list_stats_outputs[[i]] <- c(anova_outputs,lm_cpg_outputs,lm_Exp_outputs)
    names(mod_list_stats_outputs[[i]]) <- labels
    
    ### Figures 
    # Create EigenGene DNA Methylation vs. TreatmentxTime boxplot
    mod_P_A[[i]] <- ggplot(mod_list[[i]],aes(y=cpg_mean*100,x=SFV,fill=SFV)) +
      geom_boxplot() + 
      theme_cowplot() + 
      scale_fill_manual(values = col_perm[c(1,3,2,4)],
                        labels=c('Ambient\n     D9','OA 2800\n     D9',
                                 'Ambient\n     D80', 'OA 2800\n     D80')) + 
      scale_x_discrete(labels= rep("",times=4)) +
      labs(x="",y="DNA Methylation (%)",title = "") +
      theme(plot.title = element_text(hjust = 0.5),
            legend.title = element_blank(),legend.position = "none",
            axis.ticks.x = element_blank(),axis.title = element_blank())
    # Create DNA Methylation vs. EigenGene Expression
    line_val <- ifelse(mod_list_lm_cpg_expr_out[[i]]$coefficients[2,4]<0.05,1,2)
    mod_P_B[[i]] <- ggplot(mod_list[[i]],aes(y=eigen,x=cpg_mean*100,shape=SFV,colour=SFV)) + 
      geom_point(size=3) +
      geom_abline(slope = mod_list_lm_cpg_expr_out[[i]]$coefficients[2],linetype=line_val,
                  intercept = mod_list_lm_cpg_expr_out[[i]]$coefficients[1]) +
      scale_colour_manual(values = col_perm[c(1,3,2,4)],
                          labels=c('Ambient\n     D9','OA 2800\n     D9',
                                   'Ambient\n     D80', 'OA 2800\n     D80')) +
      scale_shape_manual(values = c(15,17,15,17),
                         labels=c('Ambient\n     D9','OA 2800\n     D9',
                                  'Ambient\n     D80', 'OA 2800\n     D80')) +
      theme_cowplot() + 
      labs(y="EigenGene Expression",x="DNA Methylation (%)",title="") +
      theme(plot.title = element_text(hjust = 0.5),
            legend.title = element_blank(),
            legend.position = "none",
            axis.title = element_blank())
    # Create EigenGene Expression vs. Diff. pH boxplot
    line_val <- ifelse(mod_list_lm_expr_diffpH_out[[i]]$coefficients[2,4]<0.05,1,2)
    mod_P_C[[i]] <- ggplot(mod_list[[i]],aes(y=diff_pH,x=eigen,shape=SFV,colour=SFV)) + 
      geom_point(size=3) +
      geom_abline(slope = mod_list_lm_expr_diffpH_out[[i]]$coefficients[2],
                  intercept = mod_list_lm_expr_diffpH_out[[i]]$coefficients[1],
                  linetype=line_val) +
      scale_colour_manual(values = col_perm[c(1,3,2,4)],
                          labels=c('Ambient\n     D9','OA 2800\n     D9',
                                   'Ambient\n     D80', 'OA 2800\n     D80')) +
      scale_shape_manual(values = c(15,17,15,17),
                         labels=c('Ambient\n     D9','OA 2800\n     D9',
                                  'Ambient\n     D80', 'OA 2800\n     D80')) +
      theme_cowplot() + 
      labs(y="Differential pH",x="EigenGene Expression",title="") +
      theme(plot.title = element_text(hjust = 0.5),
            legend.title = element_blank(),
            legend.position = "none",
            axis.title = element_blank())
    # Color bar
    bars <- ggplot(mod_list[[i]],aes(y=cpg_mean,x=SFV,fill=SFV)) +
      geom_boxplot() + 
      theme_cowplot() + 
      scale_fill_manual(values = col_perm[c(1,3,2,4)],
                        labels=c('Ambient\n     D9','OA 2800\n     D9','Ambient\n     D80', 'OA 2800\n     D80')) + 
      scale_x_discrete(labels= rep("",times=4)) +
      labs(x="",y="DNA Methylation (%)",title = "",
           fill="")
    legend_bars <- get_legend(bars)
    # Legend
    leg <- ggplot(mod_list[[i]],aes(y=diff_pH,x=eigen,shape=SFV,colour=SFV)) + 
      geom_point(size=3) +
      geom_abline(slope = mod_list_lm_expr_diffpH_out[[i]]$coefficients[2],
                  intercept = mod_list_lm_expr_diffpH_out[[i]]$coefficients[1]) +
      scale_colour_manual(values = col_perm[c(1,3,2,4)],
                          labels=c('Ambient\n     D9','OA 2800\n     D9','Ambient\n     D80', 'OA 2800\n     D80')) +
      scale_shape_manual(values = c(15,17,15,17),
                         labels=c('Ambient\n     D9','OA 2800\n     D9','Ambient\n     D80', 'OA 2800\n     D80')) +
      theme_cowplot() + 
      labs(y="Differential pH",x="EigenGene Expression",title="",
           shape="",colour="")
    legend_points <- get_legend(leg)
  }
  
  ### Single summary table
  mod_final <- list(Env_CpG = mod_P_A, # Plot 1 Environment vs DNA Methylation
                    CpG_Expr = mod_P_B,  # Plot 2 DNA Methylation vs Expression
                    Expr_diffpH = mod_P_C, # Plot 3 Expression vs. Diff. pH
                    ModuleSummary = mod_list, # Basic Summary Information for the module
                    annot = mod_list_full, # Table of annotated genes for each module
                    aov_env_cpg_out = mod_list_aov_env_cpg_out, # ANOVA ouput for Cpg vs TimexTreatment
                    tukey_env_cpg_out = mod_list_tukey_env_cpg_out, # Tukey output (if interaction) for ANOVA
                    lm_cpg_expr_out = mod_list_lm_cpg_expr_out, # Linear model for Module eigenExpression ~ DNA methylation
                    lm_expr_diffpH_out = mod_list_lm_expr_diffpH_out, # Linear model for Diff pH ~ Modeul eigenExpression
                    stat_summary = mod_list_stats_outputs) # Summary of stats for each analysis
                    #stat_pValue = p_list, # Pvalues for EigenGene Expression and variable associations
                    #stat_corr = corr_list) # Correlations for Eigen Expression and variable associations
  # Add module names
  for(i in 1:length(mod_final)){names(mod_final[[i]]) <- substring(targetMods,3,)}
  return(mod_final)
}

#### Running Summary Function ####

# Top three modules 
topMods <- moduleSummary(m_final,topModuleNames,modList,cpgs)
saveRDS(topMods,"results/WGCNA/WGCNA_topModSummary.RData")
allMods <- moduleSummary(m_final,rownames(modTrait_P),modList,cpgs)
saveRDS(allMods,"results/WGCNA/WGCNA_allModSummary.RData")
# Might take some time to calculate

#### Single Table Summary ####

### Table Summary Function ###
moduleTableSummary <- function(modSum) {
  # Number of genes in each module
  geneNumber <- unlist(lapply(modSum$annot,nrow))
  stat_output <- data.frame(matrix(unlist(modSum$stat_summary),ncol=length(modSum$stat_summary[[1]]),byrow = TRUE))
  names(stat_output) <- names(modSum$stat_summary[[1]])
  
  # Calculating mean expression and methylation per module
  exp_mean <- NULL
  cpg_mean <- NULL
  prop_mean <- NULL
  for (i in 1:length(modSum$ModuleSummary)){
    exp_mean <- c(exp_mean,mean(modSum$ModuleSummary[[i]]$eigen))
    cpg_mean <- c(cpg_mean,mean(modSum$ModuleSummary[[i]]$cpg_mean))
    prop_mean <- c(prop_mean,mean(modSum$ModuleSummary[[i]]$cpg_prop))
  }
  module_means <- data.frame(meanEigenExpression=exp_mean,
                             meanMethylation=cpg_mean,
                             propCpGsCovered=prop_mean)
  
  final_stats_sum <- data.frame(Module=names(modSum$Env_CpG),
                                Genes=geneNumber,
                                module_means,
                                stat_output[,c(4:6,8,7,10,12,11,14)])
}

### Creating table summary of modules

topMod_tableSummary <- moduleTableSummary(topMods)
saveRDS(topMod_tableSummary,"results/WGCNA/WGCNA_topModTableSummary.RData")
write.csv(topMod_tableSummary,"results/WGCNA/WGCNA_topModTableSummary.csv",row.names = FALSE)
allMod_tableSummary <- moduleTableSummary(allMods)
saveRDS(allMod_tableSummary,"results/WGCNA/WGCNA_allModTableSummary.RData")
write.csv(allMod_tableSummary,"results/WGCNA/WGCNA_allModTableSummary.csv",row.names = FALSE)

### READIN DATA ####
setwd("/home/downeyam/Github/AE17_Cvirginica_MolecularResponse/")
data <- allMod_tableSummary
data <- readRDS("results/WGCNA/WGCNA_allModTableSummary.RData")
topMods <-  readRDS("results/WGCNA/WGCNA_topModSummary.RData")

out2<-topMods$annot
ex<-out2$cyan
# Annotation Files for top candidates
for(i in 1:length(topMods$annot)){
  write.csv(topMods$annot[i],paste0("results/WGCNA/WGCNA_module_",names(topMods$annot)[i],"_annotation.csv"),row.names = FALSE)
}

#### Figures ####
## Labels
y.grob1 <- textGrob("Eigengene Expression",
                    gp=gpar(col="black", fontsize=15), rot=90)
x.grob1 <- textGrob("Eigengene Expression",
                    gp=gpar(col="black", fontsize=15))
y.grob2 <- textGrob("DNA Methylation (%)",
                    gp=gpar(col="black", fontsize=15), rot=90)
x.grob2 <- textGrob("DNA Methylation (%)",
                    gp=gpar(col="black", fontsize=15))
y.grob3 <- textGrob(expression(paste(Delta," pH (Total)")),
                    gp=gpar(col="black", fontsize=15), rot=90)
x.grob3 <- textGrob(expression(paste(Delta," pH (Total)")),
                    gp=gpar(col="black", fontsize=15))

## Add titles
mod_top_titles <- list()

alt_title <- c(substring(as.character(unique(topMods$ModuleSummary[[1]]$mod)),3),
               substring(as.character(unique(topMods$ModuleSummary[[2]]$mod)),3),
               substring(as.character(unique(topMods$ModuleSummary[[3]]$mod)),3))

for(i in 1:length(topMods$Env_CpG)){
  temp.title <- textGrob(substring(as.character(unique(topMods$ModuleSummary[[i]]$mod)),3),
  #temp.title <- textGrob(alt_title[i],
                          gp=gpar( col="black", fontsize=15,fontface="bold"))
  mod_top_titles[[i]] <- grid.arrange(arrangeGrob(topMods$Expr_diffpH[[i]],
                                             top=temp.title))
}

## Module Panels
mod_top <- plot_grid(plotlist = mod_top_titles,nrow=1)
mod_top_labels <- grid.arrange(arrangeGrob(mod_top,
                                           left=y.grob3,
                                           bottom=x.grob1))

mod_middle <- plot_grid(plotlist = topMods$CpG_Expr,nrow=1)
mod_middle_labels <- grid.arrange(arrangeGrob(mod_middle,
                                              left=y.grob1,bottom=x.grob2))

mod_bottom <- plot_grid(plotlist = topMods$Env_CpG,nrow=1)
mod_bottom_labels <- grid.arrange(arrangeGrob(mod_bottom,
                                              left=y.grob2))

top_legend <-  plot_grid(mod_top_labels,topMods$bars_legend,rel_widths = c(5,1))
middle_bottom <- plot_grid(mod_middle_labels,mod_bottom_labels,nrow=2)
middle_bottom_legend <- plot_grid(middle_bottom,topMods$points_legend,rel_widths = c(5,1))

mod_all <- plot_grid(top_legend,middle_bottom_legend,nrow=2,rel_heights = c(1,2))
fig7_final <-  plot_grid(NULL,mod_all,labels=c("A","B"),rel_widths = c(4,7))
fig7_final

#Alt No heatmap panel
fig7_final_alt <- plot_grid(mod_top_labels,mod_middle_labels,mod_bottom_labels,nrow=3,
                            labels=c("A","B","C"))
fig7_final_alt 
