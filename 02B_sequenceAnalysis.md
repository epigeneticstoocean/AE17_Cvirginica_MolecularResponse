# Sequence Data Analysis Workflow

## Table of Contents

1. [Data](#data)
3. [Analysis 1 - PCA and discriminant analysis of principal components](#one)
4. [Analysis 2 - Differential response](#two)
5. [Analysis 3 - Functional enrichment (GOMWU)](#three)
6. [Analysis 4 - Gene-level characterization of gene expression and DNA methylation](#four)

---

## Data <a name="data"></a>

* [**Link to scripts**](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/tree/master/src/)  
* [**Link to data**](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/tree/master/data/)  
* Reference genome: from NCBI ([GCA_002022765.4 C_virginica-3.0](https://www.ncbi.nlm.nih.gov/genome/?term=crassostrea+virginica))  

## Analysis 1 - PCA and discriminant analysis of principal components <a name="one"></a>  

Performed PCA and PERMANOVA to examine treatment and time effects on global gene expression and DNA methylation.

**Description**
We used a principal components analysis (PCA) to visualize differences among treatment and time in both global gene expression and gene body DNA methylation patterns. The gene expression PCA was based on the expression level of all genes remaining after filtering (in columns) for each individual (in rows). Gene expression level was calculated as counts per million using the cpm function from the R package edgeR (Robinson et al., 2010) and log2 transformed (i.e. log2-cpm). The DNA methylation PCA was based on all CpG loci located within gene bodies with at least 5x total coverage (in columns) for each individual (in rows).  
We used a PERMANOVA to test the null hypothesis of no effect of treatment, time, or their interaction on global gene expression and DNA methylation patterns. The PERMANOVA was based on the Manhattan distance using the adonis function in the R package vegan (v2.5-5; Oksanen 2015). Pairwise comparisons among each possible treatment and time combination were performed with the pairwise.adonis function in the pairwiseAdonis package (v. 0.3, Martinez 2019) using a FDR correction to account for multiple hypotheses.
To further investigate the specific effect of OA on genome-wide variation in both differential gene expression and DNA methylation we performed a discriminant analysis of principal components (DAPC) with the R package adegenet (Jombart et al., 2010). In brief, this method defines a discriminant function that maximally differentiates between two (or more) categories based on multidimensional sample data (e.g. expression across many genes). Here, we performed a DAPC that created a discriminant function which maximized differences between treatments in day 9 samples using the function dapc, then we predicted where samples from day 80 should fall along the discriminant function using the predict function to determine if genome-wide variation due to OA was maintained through time. 

* [script](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/src/Analyses/AE17_multivariateAnalysis.R) 

![](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/results/figures/Figure3.jpg)

## Analysis 2 - Differential response <a name="two"></a>

Performed differential expression on DNA methylation analysis on individual genes (or CpG loci within genes) to examine treatment or time effects.

**Description**
Identification of differentially expressed genes 
Differential expression among treatments for gene and transcript level expression was evaluated using a generalized linear model approach implemented in the R package limma (Ritchie et al., 2015) using treatment, time, and their interaction as fixed effects. Expression data was first TMM normalized using the calcNormFactors function and transformed into log2 counts per million (log2-cpm) using the voomWithQualityWeights function (Smyth et al., 2005). Finally, the geneDuplication function was used to account for tank as a potential experimental block effect (Oshlack et al., 2007). Site was not considered in this analysis given that it did not have a significant effect on either the phenotypic or genome-wide responses. P-values were then corrected for multiple testing using FDRtools (Strimmer, 2008). Genes with FDR ≤ 0.05 and absolute value of log2 fold ≥2 were considered differentially expressed.
Identification of differentially methylated loci 
Differential methylation was evaluated using a Bayesian binomial mixed model that included treatment, time, and their interaction as fixed factors and tank nested in shelf as a random effect using the BRMS package (Bürkner, 2017). This approach was used over more common differential methylation methods (e.g. methylkit or MACAU) because of its ability to model random effects in a generalized model framework and for its ability to improve model convergence over other generalized linear model frameworks (e.g. lme4 and MCMCglmm) using a Bayesian strategy for estimating model parameters.
In brief, we used Bayesian Regression Models (BRMs) with No-U-Turn-Sampling (NUTS) on the posterior distributions of the parameters (Bürkner, 2017). This algorithm converges much more quickly for high-dimensional models regardless of whether the priors are conjugate or not (Hoffman and Gelman 2014). To ensure that transitions after warmup did not exceed the maximum tree depth, and that there was sufficient effective sample sizes for the bulk and tail we modified several arguments including, increasing adapt_delta to 0.99, max_treedepth to 15, and number of iterations to 5000. We modeled proportion methylation as a binomial response variable with treatment level (treatment and time point) as a fixed/population effect and tank as a random/group effect with default uninformed priors. We then calculated the Bayes Factor to test three hypotheses: (i) main effect of treatment (Ambient = OA), (ii) main effect of time (day 9 = day 80), and (iii) an interaction between time and treatment. Bayes Factors were calculated with the hypothesis function in the brms package. Note that in the Bayesian framework, P-values are not calculated and therefore there is no analog to a false discovery rate correction. Hypothesis tests with log10(Bayes Factors) > 2 and 95 % credibility intervals  that did not overlap with 0 (corresponding to 20,000 posterior samples) were taken to be decisive evidence of an effect following (Kass and Raftery, 1995). For loci that had a significant interaction between time and treatment, we conducted an additional set of post-hoc hypothesis tests of individual treatment levels using the same criteria as above. 

*[Gene Expression : Differential Analysis R script](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/src/Analyses/AE17_diffMethylation.R)  
*[DNA Methlylation : Shell script for running]()  
*[DNA Methlylation : R script for performing bayesian binomial regression model]()  
*[DNA Methlylation : R script for interpreting bayesian binomial regression outputs]()  
  
## Analysis 3 - Functional enrichment (GOMWU) <a name="three"></a>  

Used the [GO-MWU](https://github.com/z0on/GO_MWU) to examine if specific genes were enriched in up or down regulated gene expression or gene level DNA methylation by treatment (at each time point).

**Description**
A functional enrichment test was conducted with GO-MWU, a rank-based gene enrichment method developed by (Wright et al., 2015), to identify gene ontology (GO) categories enriched among genes that are differentially regulated or methylated between treatments at each timepoint. We performed this analysis separately for each timepoint using the log fold change in gene expression and the difference in mean methylation among treatments. Mean methylation was calculated as the mean among all CpG loci within a gene across all individuals within a particular treatment and timepoint. Only genes with at least 5 CpG loci were considered for the analysis. Importantly, GO-MWU can handle a variety of differentiation metrics (e.g., logfold change in expression) and considers all genes, not just those that are significantly differentially expressed or methylated. This enables it to potentially detect GO categories enriched with responsive genes even when there is limited evidence of individually differentially expressed or methylated genes. 
  
Inputs for the GO-MWU analysis, include two gene list tables created using the Genebanks IDs from the C. virginica genome available on NCBI (Accession: GCA_002022765.4) along with a measure of difference (i.e. log fold change or methylation difference) and a GO terms table containing a list of all available Genebank IDs and their associated gene ontology (GO) terms. Details on how the latter file was created are outlined by Johnson et al (2019). The GO-MWU analysis was performed for both gene expression and methylation tables separately using the goStats function in R with default settings and using the gene ontology database provided by the GO-MWU repository. The analysis first clusters highly similar GO categories by combining categories that shared at least 75% of the same genes. After clustering a Mann-Whitney U test was performed to identify GO categories that were enriched with either up regulated or down regulated (or hyper- or hypo-methylated) genes. This analysis was run separately for GO categories associated with molecular function, biological process, and cellular components. A 10% FDR correction (GO-MWU default) was used to adjust for multiple comparisons. 

* [GOMWU-Folder](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/tree/master/src/Analyses/gomwu)  

## Analysis 4 - Gene-level characterization of gene expression and DNA methylation <a name="four"></a>

PCA based approach for examining correlations between gene expression and DNA mehtylation.

**Description**
To examine the relationship between gene expression, DNA methylation and gene attributes (e.g. gene length and number of exons) we used the PCA-based approach outlined in Gavery et al (2013). In the input matrix for the PCA we considered the following variables: mean gene expression over all treatments (mean log2-cpm values), the coefficient of variation (CV) in gene expression among treatment means, mean DNA methylation level over all treatments, the DNA methylation CV among treatment means, gene length, exon number per gene, and the number of CpGs per gene (script for generating each variable within the matrix are available on the github repository). Gene attributes were normalized by log transformation. Monte-Carlo randomization tests were used to determine the significance of each variable on each principal component.

* [Script](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/src/Analyses/AE17_geneAttPCA.R)

![](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/results/figures/Figure4.jpg)

