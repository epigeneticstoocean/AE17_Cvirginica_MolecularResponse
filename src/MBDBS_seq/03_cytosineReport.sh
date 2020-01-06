#!/bin/bash

dir="/shared_lab/20180226_RNAseq_2017OAExp/DNAm/trimmedSamples/singleSample_trimScript_20190719/bismark/bowtie2"
genome="/shared_lab/20180226_RNAseq_2017OAExp/DNAm/reference/genome"
out="/shared_lab/20180226_RNAseq_2017OAExp/DNAm/trimmedSamples/singleSample_trimScript_20190719"

#cd $output
counter=1

echo "Start sample # (1-24)"
read lower
echo "End sample # (1-24)"
read upper

for i in $( ls $dir/*DNAm_R1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz); do
        file=$( echo $i | rev | cut -d'/' -f 1 | rev | cut -d'_' -f 1) 
	
        if [ "$counter" -ge "$lower" ]
        then

            	if [ "$counter" -le "$upper" ]
                then

                    	echo "Sample $counter"
                        echo "File Path : $i"
                        echo "File Name : $file"

		        coverage2cytosine ${i} \
        		--genome_folder $genome \
        		--dir $out \
        		-o ${file}_CytoSummary

                fi
        fi

        counter=$((counter+1))
done

