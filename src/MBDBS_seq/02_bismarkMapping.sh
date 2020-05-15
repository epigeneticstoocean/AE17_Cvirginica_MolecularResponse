#!/bin/bash

trimmed="/shared_lab/20180226_RNAseq_2017OAExp/DNAm/trimmedSamples/singleSample_trimScript_20190719/20190719_fastqc_trim_10bp_Cvirginica_MBD"
output="/shared_lab/20180226_RNAseq_2017OAExp/DNAm/trimmedSamples/singleSample_trimScript_20190719/bismark/bowtie2"
# genome="/shared_lab/20180226_RNAseq_2017OAExp/DNAm/reference/bsgenome/wBowtie2/"

cd $output

counter=1

echo "Start sample # (1-24)"
read lower
echo "End sample # (1-24)"
read upper

for i in $( ls $trimmed/*R1*gz); do
        file=$( echo $i | rev | cut -d'/' -f 1 | rev | cut -d'_' -f 1)

        if [ "$counter" -ge "$lower" ]
        then

            	if [ "$counter" -le "$upper" ]
                then

                    	echo "Sample $counter"
                        echo "File Path : $i"
                        echo "File Name : $file"

			# Run Bismark
			bismark --non_directional -p 2\
			--score_min L,0,-0.8 \
			/shared_lab/20180226_RNAseq_2017OAExp/DNAm/reference/bsgenome/wBowtie2/ \
			-1 $trimmed/${file}*R1*gz \
			-2 $trimmed/${file}*R2*gz \
			-o $output

			# Deduplicate removal
			deduplicate_bismark -p --bam \
			$output/${file}*bam \
			--output_dir $output 

			# Samtools reorder
			samtools sort $output/${file}*deduplicated.bam \
			-o $output/${file}dedup.sorted.bam

			# Methylation Extraction
			bismark_methylation_extractor -p --bedGraph --scaffolds --counts ${file}*deduplicated.bam --multicore 20

		fi
	fi
	
	counter=$((counter+1))
done
