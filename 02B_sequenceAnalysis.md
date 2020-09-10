# Sequence Data Analysis Workflow

## Table of Contents

[Data](#data)

**Analysis**
1. [ Genome-wide analysis](#one)  
  1A.[ DNA Methylation](#oneA)  
  1B.[ RNA Seq ](#oneB)  
2. [ Differential molecular response analysis](#two)  
  2A. [ Identification of differentially methylated loci](#twoA)  
  2B. [ Identification of differentially expressed genes](#twoB)  
3. [ Gene co-expression network analysis](#three)  
4. [ Functional enrichment analysis (GOMWU)](#four)  
5. [ Gene-level characterization of gene expression and DNA methylation](#five)  

---

## Data <a name="data"></a>

[**Link to data**](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/tree/master/data/)  
  
Reference genome: from NCBI ([GCA_002022765.4 C_virginica-3.0](https://www.ncbi.nlm.nih.gov/genome/?term=crassostrea+virginica))  

## Genome-wide analysis <a name="one"></a>

### Description
The gene expression in the PCA was based on the post-filtered expression level of genes (in columns) for each individual (in rows). Gene expression level was calculated as counts per million using the cpm function from the R package `edgeR` (Robinson et al., 2010) and log2 transformed (i.e., log2-cpm). The DNA methylation PCA was based on all CpG loci located within gene bodies with at least 5x total coverage (in columns) for each individual (in rows). A PERMANOVA to test the null hypothesis of no effect of treatment, time, or their interaction on global gene expression and DNA methylation patterns. The PERMANOVA was based on the Manhattan distance using the adonis function in the R package `vegan` (v2.5-5; Dixon, 2003)

The discriminant analysis of principal components (DAPC) was performed with the R package `adegenet` (Jombart et al., 2010), using the same transcriptomic and methylomic datasets described above for the PCA. The function dapc was used to generate a discriminant function that maximized differences between treatments on day 9 samples, then we predicted where samples from day 80 should fall along the discriminant function using the predict function to determine if genome-wide variation due to OA was maintained through time. 

Lastly, we further examined the hypothesis that OA-induced changes in global DNA methylation by using a linear model to look at the effect of treatment, time, and feature (i.e., exon, intron, and intergenic region) on global methylation. Global methylation in this case was summarized as the median methylation across all CpGs within a feature for each sample.  

### DNA methylation (Figure 3) <a name="oneA"></a>

* [script](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/src/Analyses/AE17_fig3_DNAm.R) 

![](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/results/manuscript/figures/Figure3/Figure3.png)

### RNA-Seq (Figure 4) <a name="oneB"></a>

* [script](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/src/Analyses/AE17_fig4_geneExpression.R) 

![](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/results/manuscript/figures/Figure4/Figure4.png)

## Differential Expression and Methylation <a name="two"></a>

### Description 

Differential gene expression amongst treatments was evaluated using a generalized linear model approach implemented in the R package `limma` (Ritchie et al., 2015) using treatment, time, and their interaction as fixed effects. Expression data was first `'TMM'` normalized using the `calcNormFactors` function and transformed into log2 counts per million (log2-cpm) using the `voomWithQualityWeights` function (Smyth et al., 2005). Finally, the geneDuplication function was used to account for tank as a potential experimental block effect (Oshlack et al., 2007). Site was not considered in this analysis given that it did not have a significant effect on either the phenotypic or genome-wide responses. Genes with FDR ≤ 0.05 and absolute value of log2 fold ≥2 were considered differentially expressed.  

Differentially methylated loci (DML) were identified using the R package `methylKit` (Akalin et al., 2012). Only CpGs with coverage ≥ 5 for all samples were considered. Differential methylation was performed using a logistic regression approach implemented in methylKit with the functions `calculateDiffMeth` with the overdispersion argument set to `"MN"` and the default ‘slim’ method to correct p-values and the function `getMethylDiff` with a differential methylation threshold set to 50% and a q-value threshold set to 0.01.

**Differential DNA Methylation** 
* [script (towards bottom)](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/src/MBDBS_seq/04_methylationMatrix_diffMethylation_methylKit.R)

[Link to table of DMLs](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/results/manuscript/Supp/Supplemental_TableS3.3__DMLlist.csv)

**Differential Gene Expression** 

* [script](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/src/Analyses/AE17_diffExpression.R)

## Association plots (Figure 5 and 6) <a name="three"></a>

### Description

A PCA-based approach (described in Gavery and Roberts, 2013) was used to examine the relationship between gene expression, DNA methylation, and gene attributes (e.g., gene length and number of exons). The following variables were used in the input matrix for the PCA: mean gene expression over all treatments (mean log2-cpm values); the coefficient of variation (CV) in gene expression among treatment means; mean DNA methylation level over all treatments; the DNA methylation CV among treatment means, gene length, exon number per gene; and the number of CpGs per gene (script for generating each variable within the matrix are available on the github repository). Gene attributes were normalized by log transformation.  

## Gene co-expression network analysis (Figure 7) <a name="four"></a>  

### Description

A weigheted gene co-expression network analysis was performed to identify genes that exhibit similar expression patterns among individual oysters using the R package `WGCNA` (Langfelder and Horvath, 2008). First, a gene dissimilarity matrix was generated based on the log2-cpm gene expression data using first the adjacency function followed by the TOMsimilarity function in `WGCNA`. This step estimates the level of dissimilarity between each gene by considering expression across all individuals. Next, genes were hierarchically clustered based on dissimilarity using the function hclust and the `‘Ward.D2’` method for clustering (Murtagh and Legendre, 2014). Modules were determined using the cutreeDynamic function with a minimum gene membership threshold of 30. An eigenvalue for module expression (i.e., the first principle component value for each individual) was calculated for each module using moduleEigengenes. Lastly, linear regression was used to determine the association between the expression of each module (i.e., the eigenvalue of gene expression) and either mean gene methylation (calculated as the mean methylation of all CpGs among all genes within a module) or EPF response (i.e., ΔpH).


* [script]((https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/src/Analyses/AE17_fig7_WGCNAmultiComp.R)

![](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/results/manuscript/figures/Figure7/Figure7.png)

## [Functional enrichment analysis](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/tree/master/src/Analyses/gomwu) <a name="four"></a>  

A functional enrichment test was conducted with GO-MWU, a rank-based gene enrichment method developed by (Wright et al., 2015), to identify gene ontology (GO) categories enriched among genes that are differentially regulated or methylated between treatments at each time point. We performed this analysis separately for each time point using the log2-fold change in gene expression and the difference in mean methylation among treatments. Mean methylation was calculated as the mean among all CpG loci within a gene across all individuals within a particular treatment and time point. Only genes with at least 5 CpG loci were considered for the analysis to ensure mean methylation estimates were based on genes where we had at least moderate CpGs coverage. Importantly, GO-MWU can handle a variety of differentiation metrics (e.g., log2-fold change in expression) and considers all genes, not just those that are significantly differentially expressed or methylated. This enables detection of GO categories enriched with responsive genes even when there is limited evidence of individually differentially expressed or methylated genes. GO-MWU scripts and the gene ontology database were downloaded from the GO-MWU github repository (https://github.com/z0on/GO_MWU).  

  Inputs for the GO-MWU analysis include two gene list tables created using the Genebanks IDs from the C. virginica genome available on NCBI (Accession: GCA_002022765.4) along with a measure of difference (i.e., log2-fold change or methylation difference) and a table of GO terms containing a list of all available Genebank IDs and their associated gene ontology (GO) terms. Details on how the latter file was created are outlined by Johnson et al (2019) and can be found on the associated Github repository (See data availability section). The GO-MWU analysis was performed for both gene expression and methylation tables separately using the goStats function in R with default settings and using the gene ontology database provided by the GO-MWU repository. The analysis first clusters highly similar GO categories by combining categories that shared at least 75% of the same genes. After clustering, a Mann-Whitney U test was performed to identify GO categories that were enriched with either up-regulated or down-regulated (or hyper- or hypo-methylated) genes. This analysis was run separately for GO categories associated with molecular function, biological process, and cellular components. A 10% FDR correction (GO-MWU default) was used to adjust for multiple comparisons.   

## [Gene-level characterization of gene expression and DNA methylation](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/src/Analyses/AE17_fig5_DNAmvsGE.R) <a name="five"></a>  



