#!/bin/bash

## Software Versions
# STAR: STAR_2.5.3a
# RSEM: 

# Script takes trimmed RNA-seq files outputed from trimmomatic (.gz compressed format)
# and maps them to the provided genome (fasta or .fna format) and annotates with provided
# annotation file (.gtf format).

## Script uses default settings for STAR and RSEM

## WARNINGS: Make sure to adjust the number of cores used based on your computing resources.
# ALSO RSEM can take a long time. This script works well enough for only a few samples, BUT
# it is better to run RSEM in parallel if quantifying lots of samples.

### Variables and Initializations

## Number of cores
NCORE = 20

## List of samples in directory to be run
lst=(17005)
# NOTICE: This is the numbering scheme used for our samples, this might vary if you are using
# this script for other data. It ALSO assumes that there are top and bottom reads which are 
# labeled with an 'R1' and 'R2' in the file name. IF this is not the case or if the top and bottom
# reads are identified using a different schem you will need to modify the pattern on lines 
# 92,93,114,115. 
# Example of line that will need to be modified: file1=$(ls $raw/*${i}.R1.*)

# Raw sequence file directory (trimmed read files - .gz compressed)
raw="/shared_lab/20180226_RNAseq_2017OAExp/RNA/rawfiles"
# Reference genome file (.fna or fasta file)
genome="/shared_lab/20180226_RNAseq_2017OAExp/RNA/references/genome/haplotig_masked_genome.fasta"
# Annotation file (.gtf)
annotation="/shared_lab/20180226_RNAseq_2017OAExp/RNA/references/gene_annotation/KM_CV_genome_edit_Gnomon.gtf"
# Output folder (This will be created when you start the script)
output="/shared_lab/20180226_RNAseq_2017OAExp/RNA/STARTest2"

# Check if output folder already exists to prevent rewriting data
if [ -d "$output" ]; then
        echo "Directory Already Exists, please rerun with unique output directory"
        exit 1
else
        echo "Directory Created"
        mkdir "$base$output"
fi

### Create basic file structure
refs="$output/refs"
s_out="$output/STAR"
s_out_m2=$s_out"/m2"
s_out_m3=$s_out"/m3"
r_out="$output/RSEM"
mkdir $refs
mkdir $s_out
mkdir $s_out_m2
mkdir $s_out_m3
mkdir $r_out


#### Create STAR ref ##############################################
echo "Creating genome reference for STAR mapping...."
# Don't currently use this reference (I have been using the one from RSEM,
# which should be just about identitical)

# Create folder for STAR reference genome
star_ref=$refs"/STAR_ref"
mkdir $star_ref  

# Run function
STAR \
--runThreadN $NCORE \
--runMode genomeGenerate \
--genomeDir $star_ref \
--genomeFastaFiles $genome \
--sjdbGTFfile $annotation


#### RSEM from STAR ref ###########################################
echo "Creating genome reference for RSEM quantification...."

# Create folder for STAR reference genome
rsem_ref=$refs"/RSEM_ref"
mkdir $rsem_ref

# Run function
rsem-prepare-reference \
--gtf $annotation \
--star -p $NCORE \
$genome \
$rsem_ref"/RSEM"

#### STAR Mapping - first pass ####################################
echo "Starting STAR Mapping - first pass...."

# This will loop through each sample in the raw folder directory
for i in ${lst[*]};do
    echo "Star first pass .. processing sample:"$i
    outPath=$s_out_m2"/"$i

    file1=$(ls $raw/*${i}.R1.*)
    file2=$(ls $raw/*${i}.R2.*)

    # Diagnostic lines to check paths
    echo $rsem_ref
    echo $outPath

    /shared_lab/scripts/STAR \
    --runThreadN $NCORE \
    --genomeDir $rsem_ref \
    --readFilesIn $file1 $file2 \
    --outSAMtype BAM Unsorted SortedByCoordinate \
    --outFileNamePrefix $outPath \
    --readFilesCommand zcat #Used to read in compressed .gz files
done

### STAR - second pass ############################################
echo "Starting STAR Mapping - second pass...."
m2_files=$(ls $s_out_m2/*SJ.out.tab)
for i in ${lst[*]};do
    outPath=$s_out_m3"/"$i

    file1=$(ls $raw/*${i}.R1.*)
    file2=$(ls $raw/*${i}.R2.*)

    echo "Star second pass .. processing sample: "$i

    # Diagnositic echo
    echo $outPath
    echo $m2_files

## Settings
# --quantMode : 
#               TranscriptomeSAM - outputs Alignments translated into transcript coordinates
#                                  this is REQUIRED by RSEM.
#               GeneCounts - outputs a file of read counts per gene. This should correspond 
#                            to outputs generated by HtSeq
    /shared_lab/scripts/STAR \
    --runThreadN $NCORE \
    --genomeDir $rsem_ref \
    --readFilesIn $file1 $file2 \
    --outSAMtype BAM Unsorted SortedByCoordinate \
    --quantMode TranscriptomeSAM GeneCounts \
    --outFileNamePrefix $outPath \
    --readFilesCommand zcat \
    --sjdbFileChrStartEnd $m2_files # Directory from first pass with splice junctions
done

### RSEM Quantification ###########################################
echo "Starting RSEM...."

for i in ${lst[*]};do
    echo "RSEM .. processing sample: "$i
    outPath=$r_out"/"$i
    filePath=${s_out_m3}"/"${i}
    file1=$(ls $filePath*toTranscriptome.out.bam)

    # --alignments : arguement used to specify that we are supplying aligned bam files
    # --paired-end : paired end reads
    rsem-calculate-expression \
    --alignments \
    --paired-end \
    --output-genome-bam \
    -p $NCORE \
    $file1 \
    $rsem_ref"/RSEM" \
    $outPath
done