#!/bin/bash

echo "Please put in raw file directory:"
read raw
echo "Please put in name of new folder for output"
read output

base="/shared_lab/20180226_RNAseq_2017OAExp/RNA/STAR_files/"

echo "Outputs saving to : " $base$output

if [ -d "$base$output" ]; then
    echo "Directory Already Exists, please use another name"
else
    echo "Directory Created"
    mkdir "$base$output"
fi

echo "Processing the following samples: "
echo ls $raw/*.R1.*

# This will loop through each sample in the raw folder directory
for i in $( ls $raw/*.R1.* ); do
        for j in $( ls $raw/*.R2.* ); do
                file1=$( echo $i | rev | cut -d'/' -f 1 | rev | cut -d'.' -f 1)
                file2=$( echo $j | rev | cut -d'/' -f 1 | rev | cut -d'.' -f 1)

                if [ "$file1" == "$file2" ]
                then
                        echo $i and $j
                        echo RNA"$file1"_m2
                        echo yes these files match
                        /shared_lab/scripts/STAR --runThreadN 19 \
                        --genomeDir /shared_lab/20180226_RNAseq_2017OAExp/RNA/references/star_ref2 \
                        --outFilterMatchNminOverLread 0.17 --outFilterScoreMinOverLread 0.17 \
                        --readFilesIn $i $j \
                        --outSAMmapqUnique 40 \
                        --outSAMtype BAM Unsorted SortedByCoordinate \
                        --outFileNamePrefix $base$output/"$file1"_m2_ \
                        --readFilesCommand zcat
                fi
        done
done
