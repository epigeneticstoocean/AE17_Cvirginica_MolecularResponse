#!/bin/bash  

echo "Please put in raw file directory:"
read raw
#echo "Please put in name of new folder for output"
#read output

#base="/shared_lab/20180226_RNAseq_2017OAExp/RNA/STAR_files/"

#finaloutput=$base$output/m3
#m2=$base$output/m2

#echo "Outputs saving to:" $finaloutput

#echo ls "$m2" 

#if [ -d "$finaloutput" ]; then
#    echo "Directory Already Exists, please use another name"
#    exit 1
#else
#	if [ -d "$base$output/m2"]; then
#		echo "Directory Created"
#		mkdir "$finaloutput"
#	else
#		echo "First pass output doesn't exist, make sure the correct output was specified or the 1Pass script was run first."
#		exit 1
#	fi 
#fi

for i in $( ls $raw/*.R1.* ); do
        for j in $( ls $raw/*.R2.* ); do
                file1=$( echo $i | rev | cut -d'/' -f 1 | rev | cut -d'.' -f 1)
                file2=$( echo $j | rev | cut -d'/' -f 1 | rev | cut -d'.' -f 1)
                if [ "$file1" == "$file2" ]
                then
                        m2_files=$( ls /shared_lab/20180226_RNAseq_2017OAExp/RNA/STAR_files/run_20190715/m2/*m2_SJ.out.tab)
                        echo $i and $j
                        echo RNA"$file1"_m3
                        echo yes these files match
                        echo $m2_files
                        /shared_lab/scripts/STAR --runThreadN 19 \
                        --genomeDir /shared_lab/20180226_RNAseq_2017OAExp/RNA/references/star_ref2 \
                        --readFilesIn $i $j \
                        --outSAMmapqUnique 40 \
                        --outSAMtype BAM Unsorted SortedByCoordinate \
                        --quantMode TranscriptomeSAM GeneCounts --limitSjdbInsertNsj 1500000 \
                        --outFileNamePrefix /shared_lab/20180226_RNAseq_2017OAExp/RNA/STAR_files/run_20190715/m3/"$file1"_m3_ \
                        --readFilesCommand zcat \
                        --sjdbFileChrStartEnd $m2_files
                fi
        done
done
