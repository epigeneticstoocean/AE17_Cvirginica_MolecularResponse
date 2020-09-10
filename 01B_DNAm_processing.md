# MBD-BSseq Data QC, Mapping, and Quantification Pipeline

### Overview

Below are the steps to go from the raw sequence data (available on NCBI), to various methylation quantification file outputs including, methylation matrixes, sample bed files, and .RData files with and without coverage or feature filtering.

This pipeline takes advantage of a genome mapper Bismark, which is capable of alignment DNA sequence data that has been bisulfite treated. Importantly, Bismark wraps around Bowtie2 which does the actual mapping, and provides downstream commands to facilitate removing sequence duplication,and quantify cytosine and thymine coverage at different cytosine motifs across the genome (we focus on CpGs for this study).

## Table of Contents

1. [Data](#data)
2. [Brief Description and Literature on Required Tools and Scripts](#description)
3. [Step 1 - Trimming, adapter removal, and QC](#one)
4. [Step 2 - Create bisulfite treated reference genome](#two)
4. [Step 3 - Mapping with Bismark and Bowtie2](#three)
5. [Step 4 - Raw Matrices](#four)
6. [Step 5 - Filtering and Summaries](#five)

---

## Data <a name="data"></a>

* [**Link to scripts**](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/tree/master/src/MBDBS_seq)  
* [**Link to data**](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/tree/master/data)
* Reference genome: from NCBI ([GCA_002022765.4 C_virginica-3.0](https://www.ncbi.nlm.nih.gov/genome/?term=crassostrea+virginica))

## Brief Description and Literature on Required Tools and primary R packages <a name="description"></a>

**Trimming and Quality Control**

*TrimGalore!* - A flexible read trimming tool for Illumina NGS data.

* [Github](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)

**Mapping**

*Bismark* - Mapping tool wrapper for analysing bisulfite-treated data.

*[Github](https://github.com/FelixKrueger/Bismark)

*Bowtie2* - The actual mapper that is used by `Bismark`.

*[Github](https://github.com/BenLangmead/bowtie2)

**General Utility**

*samtools* - A program for basic `bed` file manipulation.  

*[Github](https://github.com/samtools/samtools)

---

## Step 1 - Trimming and Adapter Removal <a name = "one"></a>

### Overview  
Remove reads with poor mapping quality and also cut off 10bp from both the 5' and 3' regions of either strand. This will remove any adapter sequence to improve downstream mapping. This trimming follows the recommendations of bismark when the library prep was done with the pico zymo kit.

### Input

* Raw reads downloaded from NCBI bioproject.

### Code  
* [Script](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/src/MBDBS_seq/01_seqQualityTrim.sh)

Bash script for quality trimming:

Core function `trim_galore`:
```
trim_galore --paired --clip_r1 10  --clip_r2 10 \
--three_prime_clip_R1 10 \
--three_prime_clip_R2 10 \
--output_dir /path/toOutputDirectory \
--fastqc_args "--outdir /path/toFastQCOutput --threads 18" \
/path/toForwardStrand \
/path/toReverseStrand \
2> /path/to/stderr.log
```
---

## Step 2 - Create bisulfite treated reference genome <a name = "two"></a>

### Overview  
Create bisulfite converted reference genome for mapping. 

### Input

* Reference genome (see RNAseq workflow for details on downloading and prepping this).

### Code

**Command Line**
```
bismark_genome_preparation --bowtie2 --genomic_composition --parallel 10 --verbose /path/toGenomeFolder > bismark_genomePrepartion_log.txt
```
---

## Step 3 - Mapping with Bismark and Bowtie2 <a name = "three"></a>

### Overview 
The trimmed reads were mapped to the bisulfite treated reference genome (created in the previous step) in order to determine the raw counts of methylated to unmethylated cytosines at each locus. This step was saved as several outputs including a sorted .bam file and a compressed .txt file with each row as a unique CpG.

### Input

* Trimmed and QCed files from `Step 1`
* Bisulfite converted genome from `Step 2`

### Code

**Core functions for single samples**

Running `bismark` to perform mapping:
```
bismark --non_directional -p 2\
--score_min L,0,-0.8 \
path/toBisulfiteTreatedRefGenomeFolder \
-1 path/toTopTrimmedFile \
-2 path/toBottomTrimmedFile \
-o $output
```

Running `deduplicate_bismark` to remove depuplicate mappings:
```
deduplicate_bismark -p --bam \
path/toBamFile \
--output_dir $output
```

Using `samtools` to sort deduplicated bam files:
```
samtools sort path/toDeduplicated.bam \
-o path/toOutputSortedDeduplicated.bam
```

Using `bismark_methylation_extractor` to extract methylation calls:
```
bismark_methylation_extractor -p --bedGraph --scaffolds --counts path/toSortedDeduplicated.bam --multicore 20
```

Creating cytosine report (separate script) with `coverage2cytosine`:
```
coverage2cytosine path/toMethylationExtractorOutput \
--genome_folder path/toGenome \
--dir path/Direcotory \
-o path/Output 
```

* [Script that performing mapping,depublication,sorting,and methylation calls for all samples](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/src/MBDBS_seq/02_bismarkMapping.sh)
* [Script for creating full cytosine reports](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/src/MBDBS_seq/03_cytosineReport.sh)

**Disclaimer** Scripts here were hard coded and paths will need to be changed to match personal directories.

## Step 4 - Raw Matrices <a name = "four"></a>

### Overview  

This step uses a custom r script and the methylKit packaged to create a matrix for methylated and unmethylated cytosines, as well as a metaData sheet for storing information about each locus.

### Input

* MApped sample outputs from `Step 3`, specifically those generated by `coverage2cytosine`.
* [Metadata](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/tree/master/data) from folder.

### Code

* [Script](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/src/MBDBS_seq/04_methylationMatrix_diffMethylation_methylKit.R)

## Step 5  - Filtering and Summaries <a name = "five"></a>

### Overview  
Use custom r script to filter complete dataset as well as annotated and summarize by feature.

### Code
  
* [Script A](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/src/MBDBS_seq/05A_CpGCountByFeature.sh)
* [Script B](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/src/MBDBS_seq/05B_CpGIntersectionByFeature.sh)

