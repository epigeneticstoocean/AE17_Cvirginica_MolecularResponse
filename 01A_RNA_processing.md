
# RNAseq Data QC, Mapping, and Quantification Pipeline

### Overview

This pipeline takes advantage of a genome mapper STAR, which performs transcript alignment by mapping to a reference genome. Importantly, STAR is suited for the de novo discovery of splice junctions, which can be leveraged for identifying novel exons and isoforms. This pipeline couples the STAR mapper with RSEM for transcript quantification. This approach attempts to probabilistically estimate transcript abundance rather than simply count the reads.This may be beneficial for improving transcript count estimates, by probabilistically resolving reads which map to multiple genes (multimappers).  

## Table of Contents

1. [Data](#data)
2. [Brief Description and Literature on Required Tools and Scripts](#description)
3. [Step 1 - Trimming, adapter removal, and QC](#one)
4. [Step 2 - Creating STAR index](#two)
5. [Step 3 - Mapping with STAR](#three)
6. [Step 4 - Running RSEM](#four)
7. [Step 5 - Filtering, Creating DGEList Object, and Normalization (with limma-voom)](#five)
8. [Step 6 - Clustering gene expression data with WGNCA, and correlating phenotypic and environmental variables with gene clusters](#six)  

---

## Data <a name="data"></a>

* [**Link to scripts**](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/tree/master/src/RNA_seq)  
* [**Link to data**](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/tree/master/data/RNAseq)
* Reference genome: from NCBI ([GCA_002022765.4 C_virginica-3.0](https://www.ncbi.nlm.nih.gov/genome/?term=crassostrea+virginica))

## Brief Description and Literature on Required Tools and primary R packages <a name="description"></a>

**Trimming and Quality Control**

*dDocent (wrapper for trimming and QC steps)* -  dDocent pipeline uses Trimmomatic trimming tool to remove adapter and low quality sequences from the ends of reads. Within the dDocent code, it is specified to be paired-end (which is automatically recognized based on our file naming scheme), removes adapters based on thresholds for how well the adapter sequences align to reads (2:30:10; see Trimmomatic manual for more details), removes leading bases with phred quality score less than 20, removes trailing bases with phred quality score less than 20, scans the reads at a 5bp window and cuts when the average quality of the five bases is less than 10, and makes sure all reads are a minimum length after this cutting (greater than the shortest read/2).

* [Website](https://www.ddocent.com/)
* [Publication](https://peerj.com/articles/431/)

*Trimmomatic* - a  flexible read trimming tool for Illumina NGS data. Used by dDocent to trim raw RNAseq fragments and remove adapters.

* [Github](https://github.com/timflutre/trimmomatic)
* [Publication](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4103590/)

**File Conversion**

**gffread** - program to add in the conversion between different gene annotation file structures. Used here to convert from the `.gff` file format provided by NCBI to `.gtf` (preferred by STAR mapper).

* [Github](https://github.com/gpertea/gffread) : converts a `.gff` file format to `.gtf`

**Mapping**

*STAR* - fast RNA-seq aligner than can make to a reference genome and identify identify canonical as well as novel splice junctions. It will output mapped reads as `.sam` or `.bam` files, and with the `--quantMode` it can also create a tab delimited read count output (similar to HT-Seq). In addition, mapped reads can be ouputed as a `.bam` file with transcript coordinates. This can be used downstream by the transcript quantification program RSEM. 

* [Github](https://github.com/alexdobin/STAR)  
* [Publication](https://academic.oup.com/bioinformatics/article/29/1/15/272537)

**Transcript Quantification**

*RSEM* - Transcript quantifier, that can estimate counts at either the transcript (isoform) or gene level. It has a direct workflow with `STAR`, which enables a single line command for both mapping and quantification. Alternatively, it can take`STAR` outputs (specifically `.bam` files with transcript coordinates), and then perform the estimation.

* [Github](https://deweylab.github.io/RSEM/)
* [Publication](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-323)

**Gene filtering, standardization, and normalization**

*edgeR* - R package used for analyzing transcriptomic sequence data. Primarily using it to create a `DGEList` object type in R which will be used by `limma` package functions downstream. Also using it for the `cpm` function which converts within sample counts into a `count per million (cpm)`. 

* [Manual](https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)

*limma* - R package used for analyzing transcriptomic sequence data. Here we are using limma for `TMM` standardization approach to account for variable library sizes among samples, to transform our counts into `log2-cpm` using the `voom` function, account for random tank effects. Also used to perform differential expression analysis with planned contrasts in RNAseq_Analysis workflow.

* [Manual](https://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf)

---

## Step 1 - Trimming and Adapter Removal <a name = "one"></a>


Command Line:
```
./dDocent RNA.config
```
---

## Step 2 - Creating STAR index <a name="two"></a>

### Overview 
STAR performing the mapping in two primary stages. First, you need to create an index which the actuals reads are mapped too. We are creating this index using the available genome on NCBI, `.fna` file, and annotating it with a gene anotations file, `.gtf` format. 

**Additional Thoughts and Performance**
* This step only needs to be done once, unless the genome or gene annotations have been updated.
* Indexing should be relatively quick on a cluster (<10min)

**Inputs**
* [Reference genome: GCA_002022765.4 C_virginica-3.0](https://www.ncbi.nlm.nih.gov/genome/?term=crassostrea+virginica))
* [Gene annotations](https://drive.google.com/drive/u/0/folders/1KBAm4L5K2ZQ5dUOUfryk0BaB7fcA1RuL)

### Step 2.1 - File conversion

Command line code for converting from `.gff` to `.gtf`:
```
gffread my.gff -T -o my.gtf
```

### Step 2.2 - Create STAR index using oyster gene annotation file

* Script : [`STAR_genomeCreate.sh`](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/src/RNA_seq/02_STAR_genomeCreate.sh)

Command line code:
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

## Step 3 - Mapping with STAR <a name="three"></a>

### Overview
STAR maps trimmed reads to the index created in the previous steps. 

**Additional Thoughts and Performance**
* This will likely take a long time and require extensive RAM (>30GB), so will likely need to be done on a computing cluster.
* It would be a good idea to create a dettachable session via tmux as each sample takes ~2-3 hours to process.

**Input**: 

* Sample Reads (trimmed and QCed)
    * Stored as `.fq.gz` format
    * Forward Read Example : `P1_17005.R1.fq.gz`
    * Reverse Read Example : `P1_17005.R2.fq.gz`
* Index Folder (from previous step)

### **Step 3.1** : Start STAR mapping 1st Pass

* Full Script : [`STAR_1Pass_all.sh`](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/src/RNA_seq/03A_STAR_1Pass_all.sh)

Command Line:
```
downey-wall.a@comp5[references]# ./STAR_1Pass_all.sh 
Please put in raw file directory:
/pathway/to/trimmedRNASeqFiles
Please put in name of new folder for output
NAME_outputFile
```

### **Step 3.3** :  Move output files from 1st pass and Create `m3` folder for 2nd Pass

Command Line:
```
cd /pathway/to/output/folder
mkdir m2
mkdir m3
mv *m2_* m2 
```

### **Step 3.4** : Start STAR mapping 2nd Pass

* Full Script : [`STAR_2Pass_all.sh`](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/src/RNA_seq/03B_STAR_2Pass_all.sh)

Command Line:
```
downey-wall.a@comp5[references]# STAR_2Pass_all.sh
Please put in raw file directory:
/shared_lab/20180226_RNAseq_2017OAExp/RNA/rawfiles
```

**STAR command**

Core function STAR 1st pass:
```
/shared_lab/scripts/STAR --runThreadN 10 \
--genomeDir /path/toStarReferenceIndex \
--outFilterMatchNminOverLread 0.17 --outFilterScoreMinOverLread 0.17 \
--readFilesIn/path/toForwardStrand /path/toReverseStrand \
--outSAMmapqUnique 40 \
--outSAMtype BAM Unsorted SortedByCoordinate \
--outFileNamePrefix /path/toOutput \
--readFilesCommand zcat
```

Core function STAR 2nd pass:
```
/shared_lab/scripts/STAR --runThreadN 19 \
--genomeDir /path/toStarReferenceIndex \
--readFilesIn /path/toForwardStrand /path/toReverseStrand \
--outSAMmapqUnique 40 \
--outSAMtype BAM Unsorted SortedByCoordinate \
--quantMode TranscriptomeSAM GeneCounts --limitSjdbInsertNsj 1500000 \
--outFileNamePrefix /path/toOutput \
--readFilesCommand zcat \
--sjdbFileChrStartEnd /path/toSpliceJunctionFolder
```

## Step 4 - Running RSEM  <a name="four"></a>

**Creating Index folder for RSEM**

* Full Script : [`RSEM_createRefFromStar.sh`](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/src/RNA_seq/04A_RSEM_createRefFromStar.sh)

Core function `rsem-prepare-reference`: 
```
rsem-prepare-reference \
--gtf /path/toGeneAnnotation.gtf \
--star \
-p 8 \
/path/toRefGenome_GCF_002022765.2_C_virginica-3.0_genomic.fna \
/path/toOuput
```

**Performing RSEM Transcript Quantification**

* Full Script : [`RSEM_calcExp.sh`](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/src/RNA_seq/04B_RSEM_calcExp_All.sh)

Core function `rsem-calculate-expression`:

```
rsem-calculate-expression --star --paired-end \
--star-gzipped-read-file \
-p 20 \
/path/toForwardStrand \
/path/toReverseStrand \
/path/toRSEM_reference \
/path/toOutputFolder
```
## Step 5 - Filtering, Creating DGEList Object, and Normalization (with limma-voom) <a name="five"></a>

**Description**  

Takes raw rsem count estimation matrix and filters out genes that have low coverage (<1 cpm in at least 5 individuals in at least one trt/time combination), and performs normalization and transformation steps using `EdgeR` and `limma` packages.

* Full Script: [`05_filtering_CreatingDGEListObj.R`](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/src/RNA_seq/05_filtering_CreatingDGEListObj.R)

## Step 6 - Clustering gene expression data with WGNCA, and correlating phenotypic and environmental variables with gene clusters <a name="six"></a>

**Description**
  
A weighted gene co-expression network analysis was performed to identify genes that exhibit similar expression patterns among individual oysters using the R package WGCNA (Langfelder and Horvath, 2008). We followed a standard WGCNA pipeline for clustering, association testing, and creating WGCNA objects. This analysis was performed on the 22 individuals that remained after excluding individuals that were identified as outliers in either the gene expression (17005) and DNA methylation (17099) data. First, a gene dissimilarity matrix was generated based on the log2-cpm gene expression data using first the adjacency function followed by the TOMsimilarity function in WGCNA. This step estimates the level of dissimilarity between each gene by considering expression across all individuals. Next, genes were hierarchically clustered based on dissimilarity using the function hclust and the ‘Ward.D2’ method for clustering (Murtagh and Legendre, 2014). There objects created hear were used in the analysis step to generate figure 7.

* Full Script: [`06_CreatingWGCNAObj.R`](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/src/RNA_seq/06_CreatingWGCNAObj.R)
