#!/bin/bash

# Script intersects CpG bed tracks with the genome track file, counting the number of CpGs from
# each track in each gene. It also calculates mean summary values for several CpG related stats
# for each gene (i.e. summarizes across CpGs within each gene)

#Path variables and date
saveDate="20200130"
outputDir="/shared_lab/20180226_RNAseq_2017OAExp/DNAm/processed_samples/09_CpG_summary/"
CpGpath="/shared_lab/20180226_RNAseq_2017OAExp/DNAm/processed_samples/03_CytoSummaries/methylKit_outputs/"
genomeTrack="/shared_lab/20180226_RNAseq_2017OAExp/genome-feature-tracks/"
pathBedTools="/home/downey-wall.a/software/bedtools2/bin/"

# Different CpG bed files
all=${CpGpath}20190129_CpG_all_Destrand.bed
cov5=${CpGpath}20190129_CpG_5xCov_Destrand.bed
DML_trt=${CpGpath}20190129_DML_all_Destrand_5x_50_q01.bed # DMLs from methylKit
DML_tp9=${CpGpath}20190129_DML_tp9_Destrand_5x_50_q01.bed # DMLs from methylKit
DML_tp80=${CpGpath}20190129_DML_tp80_Destrand_5x_50_q01.bed # DMLs from methylKit
# Genome track bed file
geneL=${genomeTrack}C_virginica-3.0_Gnomon_gene_sorted_yrv.bed

outputFile=${outputDir}${saveDate}_CpGbyGeneSummary
mkdir -p ${outputFile}

#### Intersect and Count all CpGs for each gene ####
echo "Intersect and Count all CpGs for each gene... "

${pathBedTools}bedtools intersect \
-wa -wb \
-a ${geneL} \
-b ${all} \
> ${outputFile}/${saveDate}_temp_allIntersect.txt

${pathBedTools}bedtools groupby \
-i ${outputFile}/${saveDate}_temp_allIntersect.txt \
-g 1-3 \
-c 7 \
-o count \
> ${outputFile}/${saveDate}_allCpG_Count.txt

#### Intersect and Count and Summarize Methylation for all CpGs with 5x coverage for each gene ####
echo "Counting up all CpGs with 5x coverage in genes ... "

${pathBedTools}bedtools intersect \
-wa -wb \
-a ${geneL} \
-b ${cov5} \
> ${outputFile}/${saveDate}_temp_cov5Intersect.txt

${pathBedTools}bedtools groupby \
-i ${outputFile}/${saveDate}_temp_cov5Intersect.txt \
-g 1-3 \
-c 7 \
-o count \
> ${outputFile}/${saveDate}_cov5CpG_Count.txt

# Summarize the methylation for CpGs with coverage 
${pathBedTools}bedtools groupby \
-i ${outputFile}/${saveDate}_temp_cov5Intersect.txt \
-g 1-3 \
-c 8,9,10,11,12,13,14,15,16,17 \
-o mean \
> ${outputFile}/${saveDate}_cov5CpG_Summarize.txt

#### Intersect and Count all DMLs (for treatment) for each gene ####
echo "Intersect and Count all DMLs (for treatment) for each gene... "

${pathBedTools}bedtools intersect \
-wa -wb \
-a ${geneL} \
-b ${DML_trt} \
> ${outputFile}/${saveDate}_temp_DMLTrtIntersect.txt

${pathBedTools}bedtools groupby \
-i ${outputFile}/${saveDate}_temp_DMLTrtIntersect.txt \
-g 1-3 \
-c 7 \
-o count \
> ${outputFile}/${saveDate}_DMLTrt_Count.txt

#### Intersect and Count all DMLs (for treatment at tp9) for each gene ####
echo "Intersect and Count all DMLs (for treatment at tp9) for each gene... "

${pathBedTools}bedtools intersect \
-wa -wb \
-a ${geneL} \
-b ${DML_tp9} \
> ${outputFile}/${saveDate}_temp_DML_tp9Intersect.txt

${pathBedTools}bedtools groupby \
-i ${outputFile}/${saveDate}_temp_DML_tp9Intersect.txt \
-g 1-3 \
-c 7 \
-o count \
> ${outputFile}/${saveDate}_DML_tp9_Count.txt

#### Intersect and Count all DMLs (for treatment at tp80) for each gene ####
echo "Intersect and Count all DMLs (for treatment at tp80) for each gene... "

${pathBedTools}bedtools intersect \
-wa -wb \
-a ${geneL} \
-b ${DML_tp80} \
> ${outputFile}/${saveDate}_temp_DML_tp80Intersect.txt

${pathBedTools}bedtools groupby \
-i ${outputFile}/${saveDate}_temp_DML_tp80Intersect.txt \
-g 1-3 \
-c 7 \
-o count \
> ${outputFile}/${saveDate}_DML_tp80_Count.txt

# We are also going to generate a gene to exon intersect and count to count the
# number of exons in each gene.
echo "Counting all exons in each gene... "

${pathBedTools}bedtools intersect \
-wa -wb \
-a ${geneL} \
-b ${exonL} \
> ${outputFile}/${saveDate}_temp_exon_tp9Intersect.txt

${pathBedTools}bedtools groupby \
-i ${outputFile}/${saveDate}_temp_exon_tp9Intersect.txt \
-g 1-3 \
-c 7 \
-o count \
> ${outputFile}/${saveDate}_exonInGene_Count.txt

rm -rf ${outputFile}/*temp*

