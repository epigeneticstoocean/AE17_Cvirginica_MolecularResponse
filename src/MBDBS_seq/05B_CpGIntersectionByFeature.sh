#!/bin/bash

# Script takes the CpG bed files generated from methyKit work up of the Cyto_summary outputs from bismark and
# intersect them with the major genomic features (genes, exons,introns, intergenic,mRNA) using genome tracks generated
# by Yaamini. It does this for all CpGs, those with at least 5x coverage, and those that were diff. mehtylated by trt
# among both timepoints, and then for each timepoint separately.

## Setting variable paths

# Date used when auto generating your file outputs
# Output Path
path="/shared_lab/20180226_RNAseq_2017OAExp/DNAm/processed_samples/09_CpG_summary/"
# Date header on each file
saveDate="20200129"
# Path to where the bed files are located
CpGPath="/shared_lab/20180226_RNAseq_2017OAExp/DNAm/processed_samples/03_CytoSummaries/methylKit_outputs/"
# Path to the genome track files
featurePath="/shared_lab/20180226_RNAseq_2017OAExp/genome-feature-tracks/"

# Differentially methylated loci lists (as .bed files)
CpG_all="20190129_CpG_all_Destrand.bed"
CpG_cov5="20190129_CpG_5xCov_Destrand.bed"
DML_t="20190129_DML_all_Destrand_5x_50_q01.bed"
DML_9="20190129_DML_tp9_Destrand_5x_50_q01.bed"
DML_80="20190129_DML_tp80_Destrand_5x_50_q01.bed"
# Hyper vs. Hypo
DML_t_hyper="20190129_DML_all_Destrand_5x_50_q01_hyper.bed"
DML_t_hypo="20190129_DML_all_Destrand_5x_50_q01_hypo.bed"
DML_9_hyper="20190129_DML_tp9_Destrand_5x_50_q01_hyper.bed"
DML_9_hypo="20190129_DML_tp9_Destrand_5x_50_q01_hypo.bed"
DML_80_hyper="20190129_DML_tp80_Destrand_5x_50_q01_hyper.bed"
DML_80_hypo="20190129_DML_tp80_Destrand_5x_50_q01_hypo.bed"
# Feature lists (as .bed files)
geneGFF="C_virginica-3.0_Gnomon_gene_yrv.gff3"
exonL="C_virginica-3.0_Gnomon_exon_sorted_yrv.bed"
intronL="C_virginica-3.0_Gnomon_intron_yrv.bed"
geneL="C_virginica-3.0_Gnomon_gene_sorted_yrv.bed"
interL="C_virginica-3.0_Gnomon_intergenic_yrv.bed"
mRNAL="C_virginica-3.0_Gnomon_mRNA_yrv.bed"

#### Overlap of exons with genes ####
echo "Identifying ${exonL} in ${geneL} ..." 
echo "Identifying ${exonL} in ${geneL} ..." > ${path}${saveDate}_overlappingSummary.txt
# Exons in Genes
bedtools intersect \
-wa -wb \
-a ${featurePath}${geneL} \
-b ${featurePath}${exonL} \
> ${path}${saveDate}_exonInGeneSummary.txt
exonInGene="${path}${saveDate}_exonInGeneSummary.txt"

# Annotate
bedtools intersect \
-wa -wb \
-a ${exonInGene} \
-b ${featurePath}${geneGFF} \
> ${path}${saveDate}_exonInGeneSummary_annotated.txt

#### Examine how the CpGs overlap with different features ####
echo "Identifying CpGs in ${CpG_cov5} ..." 
echo "Identifying CpGs in ${CpG_cov5} ..." >> ${path}${saveDate}_overlappingSummary.txt
# Exons
bedtools intersect \
-wa -wb \
-b ${CpGPath}${CpG_cov5} \
-a ${featurePath}${exonL} \
> ${path}${saveDate}_CpG_cov5_Exon.txt
echo "that overlap with ${exonL}..." >> ${path}${saveDate}_overlappingSummary.txt

# Introns
bedtools intersect \
-wa -wb \
-b ${CpGPath}${CpG_cov5} \
-a ${featurePath}${intronL} \
> ${path}${saveDate}_CpG_cov5_Intron.txt
echo "that overlap with ${intronL}..." >> ${path}${saveDate}_overlappingSummary.txt

# Genes
bedtools intersect \
-wa -wb \
-b ${CpGPath}${CpG_cov5} \
-a ${featurePath}${geneL} \
> ${path}${saveDate}_CpG_cov5_gene.txt
echo "that overlap with ${geneL}..." >> ${path}${saveDate}_overlappingSummary.txt

# Intergenic
bedtools intersect \
-wa -wb \
-b ${CpGPath}${CpG_cov5} \
-a ${featurePath}${interL} \
> ${path}${saveDate}_CpG_cov5_Intergenic.txt
echo "that overlap with ${interL}..." >> ${path}${saveDate}_overlappingSummary.txt

# mRNA
bedtools intersect \
-wa -wb \
-b ${CpGPath}${CpG_cov5} \
-a ${featurePath}${mRNAL} \
> ${path}${saveDate}_CpG_cov5_mRNA.txt
echo "that overlap with ${mRNAL}..." >> ${path}${saveDate}_overlappingSummary.txt

## Examine how the CpGs overlap with different features
echo "Identifying all CpGs in ${CpG_all} ..." >> ${path}${saveDate}_overlappingSummary.txt
# Exons 
bedtools intersect \
-wa -wb \
-b ${CpGPath}${CpG_all} \
-a ${featurePath}${exonL} \
> ${path}${saveDate}_CpG_all_Exon.txt
echo "that overlap with ${exonL}..." >> ${path}${saveDate}_overlappingSummary.txt

# Introns
bedtools intersect \
-wa -wb \
-b ${CpGPath}${CpG_all} \
-a ${featurePath}${intronL} \
> ${path}${saveDate}_CpG_all_Intron.txt
echo "that overlap with ${intronL}..." >> ${path}${saveDate}_overlappingSummary.txt

# Genes
bedtools intersect \
-wa -wb \
-b ${CpGPath}${CpG_all} \
-a ${featurePath}${geneL} \
> ${path}${saveDate}_CpG_all_gene.txt
echo "that overlap with ${geneL}..." >> ${path}${saveDate}_overlappingSummary.txt

# Intergenic
bedtools intersect \
-wa -wb \
-b ${CpGPath}${CpG_all} \
-a ${featurePath}${interL} \
> ${path}${saveDate}_CpG_all_Intergenic.txt
echo "that overlap with ${interL}..." >> ${path}${saveDate}_overlappingSummary.txt

# mRNA
bedtools intersect \
-wa -wb \
-b ${CpGPath}${CpG_all} \
-a ${featurePath}${mRNAL} \
> ${path}${saveDate}_CpG_all_mRNA.txt
echo "that overlap with ${mRNAL}..." >> ${path}${saveDate}_overlappingSummary.txt


## Examine how the DMLs (DML_t) for treatment overlap with different features
echo "Identifying all CpGs in ${DML_t} ..." >> ${path}${saveDate}_overlappingSummary.txt
# Exons 
bedtools intersect \
-wa -wb \
-b ${CpGPath}${DML_t} \
-a ${featurePath}${exonL} \
> ${path}${saveDate}_DML_all_Exon.txt
echo "that overlap with ${exonL}..." >> ${path}${saveDate}_overlappingSummary.txt

# Introns
bedtools intersect \
-wa -wb \
-b ${CpGPath}${DML_t} \
-a ${featurePath}${intronL} \
> ${path}${saveDate}_DML_all_Intron.txt
echo "that overlap with ${intronL}..." >> ${path}${saveDate}_overlappingSummary.txt

# Genes
bedtools intersect \
-wa -wb \
-b ${CpGPath}${DML_t} \
-a ${featurePath}${geneL} \
> ${path}${saveDate}_DML_all_gene.txt
echo "that overlap with ${geneL}..." >> ${path}${saveDate}_overlappingSummary.txt

# Intergenic
bedtools intersect \
-wa -wb \
-b ${CpGPath}${DML_t} \
-a ${featurePath}${interL} \
> ${path}${saveDate}_DML_all_Intergenic.txt
echo "that overlap with ${interL}..." >> ${path}${saveDate}_overlappingSummary.txt

# mRNA
bedtools intersect \
-wa -wb \
-b ${CpGPath}${DML_t} \
-a ${featurePath}${mRNAL} \
> ${path}${saveDate}_DML_all_mRNA.txt
echo "that overlap with ${mRNAL}..." >> ${path}${saveDate}_overlappingSummary.txt

## Examine how the DMLs (DML_9) for treatment overlap with different features
echo "Identifying all CpGs in ${DML_9} ..." >> ${path}${saveDate}_overlappingSummary.txt
# Exons 
bedtools intersect \
-wa -wb \
-b ${CpGPath}${DML_9} \
-a ${featurePath}${exonL} \
> ${path}${saveDate}_DML_tp9_Exon.txt
echo "that overlap with ${exonL}..." >> ${path}${saveDate}_overlappingSummary.txt

# Introns
bedtools intersect \
-wa -wb \
-b ${CpGPath}${DML_9} \
-a ${featurePath}${intronL} \
> ${path}${saveDate}_DML_tp9_Intron.txt
echo "that overlap with ${intronL}..." >> ${path}${saveDate}_overlappingSummary.txt

# Genes
bedtools intersect \
-wa -wb \
-b ${CpGPath}${DML_9} \
-a ${featurePath}${geneL} \
> ${path}${saveDate}_DML_tp9_gene.txt
echo "that overlap with ${geneL}..." >> ${path}${saveDate}_overlappingSummary.txt

# Intergenic
bedtools intersect \
-wa -wb \
-b ${CpGPath}${DML_9} \
-a ${featurePath}${interL} \
> ${path}${saveDate}_DML_tp9_Intergenic.txt
echo "that overlap with ${interL}..." >> ${path}${saveDate}_overlappingSummary.txt

# mRNA
bedtools intersect \
-wa -wb \
-b ${CpGPath}${DML_9} \
-a ${featurePath}${mRNAL} \
> ${path}${saveDate}_DML_tp9_mRNA.txt
echo "that overlap with ${mRNAL}..." >> ${path}${saveDate}_overlappingSummary.txt

## Examine how the DMLs (DML_80) for treatment overlap with different features
echo "Identifying all CpGs in ${DML_80} ..." >> ${path}${saveDate}_overlappingSummary.txt
# Exons 
bedtools intersect \
-wa -wb \
-b ${CpGPath}${DML_80} \
-a ${featurePath}${exonL} \
> ${path}${saveDate}_DML_tp80_Exon.txt
echo "that overlap with ${exonL}..." >> ${path}${saveDate}_overlappingSummary.txt

# Introns
bedtools intersect \
-wa -wb \
-b ${CpGPath}${DML_80} \
-a ${featurePath}${intronL} \
> ${path}${saveDate}_DML_tp80_Intron.txt
echo "that overlap with ${intronL}..." >> ${path}${saveDate}_overlappingSummary.txt

# Genes
bedtools intersect \
-wa -wb \
-b ${CpGPath}${DML_80} \
-a ${featurePath}${geneL} \
> ${path}${saveDate}_DML_tp80_gene.txt
echo "that overlap with ${geneL}..." >> ${path}${saveDate}_overlappingSummary.txt

# Intergenic
bedtools intersect \
-wa -wb \
-b ${CpGPath}${DML_80} \
-a ${featurePath}${interL} \
> ${path}${saveDate}_DML_tp80_Intergenic.txt
echo "that overlap with ${interL}..." >> ${path}${saveDate}_overlappingSummary.txt

# mRNA
bedtools intersect \
-wa -wb \
-b ${CpGPath}${DML_80} \
-a ${featurePath}${mRNAL} \
> ${path}${saveDate}_DML_tp80_mRNA.txt
echo "that overlap with ${mRNAL}..." >> ${path}${saveDate}_overlappingSummary.txt
