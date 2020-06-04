#!/bin/bash

for i in $( ls /shared_lab/20180226_RNAseq_2017OAExp/RNA/rawfiles/*.R1.* ); do
        for j in $( ls /shared_lab/20180226_RNAseq_2017OAExp/RNA/rawfiles/*.R2.* ); do
                file1=$( echo $i | rev | cut -d'/' -f 1 | rev | cut -d'.' -f 1)
                file2=$( echo $j | rev | cut -d'/' -f 1 | rev | cut -d'.' -f 1)

                if [ "$file1" == "$file2" ]
                then
                        echo "Running Sample" $file1
                        rsem-calculate-expression --star --paired-end \
                        --star-gzipped-read-file \
                        -p 20 \
                        $i \
                        $j \
                        /shared_lab/20180226_RNAseq_2017OAExp/RNA/references/RSEM_refFromGFF/RSEM_ref \
                        /shared_lab/20180226_RNAseq_2017OAExp/RNA/STAR_files/run_20190707_RSEMdirectWithSTAR/RSEM_$file1_
                fi
        done
done

