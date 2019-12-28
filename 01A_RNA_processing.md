
# RNAseq Data QC, Mapping, and Quantification Pipeline

### Overview

This pipeline takes advantage of a genome mapper STAR, which performs transcript alignment by mapping to a reference genome. Importantly, STAR is suited for the de novo discovery of splice junctions, which can be leveraged for identifying novel exons and isoforms. This pipeline couples the STAR mapper with RSEM for transcript quantification. This approach attempts to probabilistically estimate transcript abundance rather than simply count the reads.This may be beneficial for improving transcript count estimates, by probabilistically resolving reads which map to multiple genes (multimappers).  

## Table of Contents

1. [Data](#zero)
2. [Brief Description and Literature on Required Tools and Scripts](#one)
3. [Step 1 - Creating STAR index](#two)
4. [Step 2 - Mapping with STAR](#three)
5. [Step 3 - Running RSEM](#four)

## Data <a name="zero"></a>

[**Link to scripts**](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/tree/master/src/RNAseq)  
[**Link to data**](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/tree/master/data/RNAseq)  

## Brief Description and Literature on Required Tools and Scripts <a name="one"></a>

**Mapping**

*STAR* - fast RNA-seq aligner than can make to a reference genome and identify identify canonical as well as novel splice junctions. It will output mapped reads as `.sam` or `.bam` files, and with the `--quantMode` it can also create a tab delimited read count output (similar to HT-Seq). In addition, mapped reads can be ouputed as a `.bam` file with transcript coordinates. This can be used downstream by the transcript quantification program RSEM. 

* [Github](https://github.com/alexdobin/STAR)  
* [Publication](https://academic.oup.com/bioinformatics/article/29/1/15/272537)

**Transcript Quantification**

*RSEM* - Transcript quantifier, that can estimate counts at either the transcript (isoform) or gene level. It has a direct workflow with `STAR`, which enables a single line command for both mapping and quantification. Alternatively, it can take`STAR` outputs (specifically `.bam` files with transcript coordinates), and then perform the estimation.

* [Github](https://deweylab.github.io/RSEM/)
* [Publication](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-323)

## Step 1 - Creating STAR index <a name="two"></a>

### Overview 
STAR performing the mapping in two primary stages. First, you need to create an index which the actuals reads are mapped too. We are creating this index using the available genome on NCBI, `.fna` file, and annotating it with a gene anotations file, `.gtf` format. 

**Additional Thoughts and Performance**
* This step only needs to be done once, unless the genome or gene annotations have been updated.
* Indexing should be relatively quick on a cluster (<10min)

**Inputs**
* Reference genome: from NCBI ([GCA_002022765.4 C_virginica-3.0](https://www.ncbi.nlm.nih.gov/genome/?term=crassostrea+virginica))
* Gene annotations: created by Kevin Johnson during the frogger workshop [LINK](https://drive.google.com/drive/u/0/folders/1KBAm4L5K2ZQ5dUOUfryk0BaB7fcA1RuL), which I converted from a `.gff` file formate to `.gtf` using [`gffread`](https://github.com/gpertea/gffread).

**Outputs**
* Reference Star Folder: `/shared_lab/20180226_RNAseq_2017OAExp/RNA/references/mapping_indexes/STAR_gnomon`

### Coding and Scripts

Command line code for converting from `.gff` to `.gtf`:
```
gffread my.gff -T -o my.gtf
```

Command line code for creating STAR index for oyster:
```
downey-wall.a@comp5[references]# STAR_genomeCreate.sh 
/shared_lab/20180226_RNAseq_2017OAExp/RNA/scripts/STAR_scripts/STAR_genomeCreate.sh: line 1: !#/bin/bash: No such file or directory
Please put in the base directory:
/shared_lab/20180226_RNAseq_2017OAExp/RNA/references/
Please put in the output folder name:
star_ref2
Outputs saving to :  /shared_lab/20180226_RNAseq_2017OAExp/RNA/references/star_ref2
Directory Created
Select genome file (.fna format, should include entire path)
/shared_lab/20180226_RNAseq_2017OAExp/RNA/references/genome/GCF_002022765.2_C_virginica-3.0_genomic.fna
Select gene annotation file (.gtf, should includ entire path)
/shared_lab/20180226_RNAseq_2017OAExp/RNA/references/gene_annotation/KM_CV_genome.gtf 
```

Bash code for bash script `STAR_genomeCreate.sh`:
```
!#/bin/bash

# Prompts user to input a path where the files created will be stored and the name of the new folder
echo "Please put in the base directory:"
read base
echo "Please put in the output folder name:"
read output

echo "Outputs saving to : " $base$output

if [ -d "$base$output" ]; then
        echo "Directory Already Exists, please rerun with unique output directory"
        exit 1
else
    	echo "Directory Created"
        mkdir "$base$output"
fi

# User selects the genome file
echo "Select genome file (.fna format, should include entire path)"
read genome

#User selects the gene annotation file (in the gtf format)
echo "Select gene annotation file (.gtf, should includ entire path)"
read gene_annotation

# Run actual indexing step
STAR --runThreadN 32 \
--runMode genomeGenerate \
--genomeDir $base$output \
--genomeFastaFiles $genome \
--sjdbGTFfile $gene_annotation
```

Output code if run successfully:
```
Jul 15 12:23:09 ..... started STAR run
Jul 15 12:23:09 ... starting to generate Genome files
Jul 15 12:23:22 ... starting to sort Suffix Array. This may take a long time...
Jul 15 12:23:25 ... sorting Suffix Array chunks and saving them to disk...
Jul 15 12:24:53 ... loading chunks from disk, packing SA...
Jul 15 12:25:29 ... finished generating suffix array
Jul 15 12:25:29 ... generating Suffix Array index
Jul 15 12:27:55 ... completed Suffix Array index
Jul 15 12:27:55 ..... processing annotations GTF
Jul 15 12:28:06 ..... inserting junctions into the genome indices
Jul 15 12:29:46 ... writing Genome to disk ...
Jul 15 12:29:48 ... writing Suffix Array to disk ...
Jul 15 12:30:07 ... writing SAindex to disk
Jul 15 12:30:12 ..... finished successfully
```

--- 

## Step 2 - Mapping with STAR <a name="three"></a>

### Overview
STAR maps trimmed reads to the index created in the previous steps. 

**Additional Thoughts and Performance**
* This will likely take a long time and require extensive RAM (>30GB), so will likely need to be done on a computing cluster.
* It would be a good idea to create a dettachable session via tmux as each sample takes ~2-3 hours to process.

**Input**: 

* Sample Reads
    * Stored as `.fq.gz` format
    * Folder Name: `/shared_lab/20180226_RNAseq_2017OAExp/RNA/rawfiles`
    * Forward Read Example : `P1_17005.R1.fq.gz`
    * Reverse Read Example : `P1_17005.R2.fq.gz`
    * This are fastq files that contains reads that have been trimmed and undergone basic quality control following the dDocent pipeline and trimmomatic.

* Index Folder (from previous step)
    * Folder path: `/shared_lab/20180226_RNAseq_2017OAExp/RNA/references/mapping_indexes/STAR_gnomon`

**Output**:
* Base Directory for STAR output files: `/shared_lab/20180226_RNAseq_2017OAExp/RNA/STAR_files/final_data`
        * First pass results : `/m2`
        * Second pass results : `/m3`

### Coding and Scripts

**Step 2.1** : Create TMUX session (Optional: do this once you've ssh'ed into the cluster not before):
* This will take a long time to run, so it will be goo

Command Line:
```
tmux new-session -s NAME_SESSION_HERE
```

---

**Step 2.2** : Start STAR mapping 1st Pass

Command Line:
```
downey-wall.a@comp5[references]# STAR_1Pass_all.sh 
Please put in raw file directory:
/shared_lab/20180226_RNAseq_2017OAExp/RNA/rawfiles
Please put in name of new folder for output
run_20190715
```
---

**Step 2.3** :  Move output files from 1st pass and Create `m3` folder for 2nd Pass
```
cd /pathway/to/output/folder
mkdir m2
mkdir m3
mv *m2_* m2 
```
* This step is need because the second pass script is currently not super flexible to user input. In the future this step  will be removed.

---

**Step 2.4** : Start STAR mapping 2nd Pass

Command Line:
```
downey-wall.a@comp5[references]# STAR_1Pass_all.sh
Please put in raw file directory:
/shared_lab/20180226_RNAseq_2017OAExp/RNA/rawfiles
```

Dettach TMUX session to prevent accidental pipe breaking, using hot keys. 
`Ctrl` + `Shift` + `B` + `D`

**BASH Scripts**

Bash code for bash script for first pass : `STAR_1Pass_all.sh`:
```
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
                        /shared_lab/scripts/STAR --runThreadN 10 \
                        --genomeDir /shared_lab/20180226_RNAseq_2017OAExp/RNA/references/RSEM_gnomon \
                        --outFilterMatchNminOverLread 0.17 --outFilterScoreMinOverLread 0.17 \
                        --readFilesIn $i $j \
                        --outSAMmapqUnique 40 \
                        --outSAMtype BAM Unsorted SortedByCoordinate \
                        --outFileNamePrefix $base$output/"$file1"_m2_ \
                        --readFilesCommand zcat
                fi
        done
done
```

Bash Code for the 2nd Pass, `STAR_2Pass_all.sh`:
```
#!/bin/bash

echo "Please put in raw file directory:"
read raw

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
```
