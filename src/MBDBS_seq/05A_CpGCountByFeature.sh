#!/bin/bash

# Script counts the occurence of CpGs within different file subsets (all,5x coverage, DMLs etc)
# for each genomic feature (gene, exon, etc.). 
# For example, it will create a file that counts the number of CpGs (from all possible CpGs or one of the subsets) 
# that exist in each gene. So you will have a list of genes with CpGs, and the number 
# of CpGs that exist within that gene.  

## Setting variable paths

# Date used when auto generating your file outputs
path="/shared_lab/20180226_RNAseq_2017OAExp/DNAm/processed_samples/09_CpG_summary/"
pathBedTools="/home/downey-wall.a/software/bedtools2/bin/"
saveDate="20200130"

cd ${path}

mkdir -p ${saveDate}_groupbySummaries
out=${saveDate}_groupbySummaries/

for file in *CpG_all*
do
echo "Grouping ${file}..."
        if wc -l ${file} > 0 ; then
                feature=`echo ${file}|cut -f 4 -d '_' |cut -f 1 -d '.'`
                ${pathBedTools}bedtools groupby \
                -i ${file} \
                -g 1-3 \
                -c 7 \
                -o count \
                > ${out}CpG_all_${feature}.txt
        fi
done

for file in *CpG_cov5*
do
echo "Grouping ${file}..."
        if wc -l ${file} > 0 ; then
                feature=`echo ${file}|cut -f 4 -d '_' |cut -f 1 -d '.'`
                ${pathBedTools}bedtools groupby \
                -i ${file} \
                -g 1-3 \
                -c 7 \
                -o count \
                > ${out}CpG_cov5_${feature}.txt
        fi
done


for file in *DML_all*
do
echo "Grouping ${file}..."
        if wc -l ${file} > 0 ; then
                feature=`echo ${file}|cut -f 4 -d '_' |cut -f 1 -d '.'`
                ${pathBedTools}bedtools groupby \
                -i ${file} \
                -g 1-3 \
                -c 7 \
                -o count \
                > ${out}DML_all_${feature}.txt
        fi
done


for file in *DML_tp9*
do
echo "Grouping ${file}..."
        if wc -l ${file} > 0 ; then
                feature=`echo ${file}|cut -f 4 -d '_' |cut -f 1 -d '.'`
                ${pathBedTools}bedtools groupby \
                -i ${file} \
                -g 1-3 \
                -c 7 \
                -o count \
                > ${out}DML_tp9_${feature}.txt
        fi
done


for file in *DML_tp80*
do
echo "Grouping ${file}..."
        if wc -l ${file} > 0 ; then
                feature=`echo ${file}|cut -f 4 -d '_' |cut -f 1 -d '.'`
                ${pathBedTools}bedtools groupby \
                -i ${file} \
                -g 1-3 \
                -c 7 \
                -o count \
                > ${out}DML_tp80_${feature}.txt
        fi
done

# This additional code will take all the count file for each subset of CpGs (all,cov5,etc) and
# and combine them into a single file per genomic feature.
 
echo "Summarizing all CpG subsets by feature..."

cd ${path}${out}

ftr=("gene" "Exon" "Intron" "Intragenic")
for i in ${ftr[@]}
do
  t1="CpG_all_${i}.txt"
  t2="CpG_cov5_${i}.txt"
  t3="DML_all_${i}.txt"
  t4="DML_tp9_${i}.txt"
  t5="DML_tp80_${i}.txt"

  bedtools intersect -wa -wb \
  -loj -a ${t1} -b ${t2} > temp_${i}_1.txt

  bedtools intersect -wa -wb \
  -loj -a temp_${i}_1.txt -b ${t3} > temp_${i}_2.txt

  bedtools intersect -wa -wb \
  -loj -a temp_${i}_2.txt -b ${t4} > temp_${i}_3.txt

  bedtools intersect -wa -wb \
  -loj -a temp_${i}_3.txt -b ${t5} > fullFeatureCountSummary_${i}.txt
done

rm -rf temp*

