# MBD-BSseq Data QC, Mapping, and Quantification Pipeline

### Overview

Below are the steps to go from the raw sequence data (available on NCBI), to various methylation quantification file outputs including, methylation matrixes, sample bed files, and .RData files with and without coverage or feature filtering.

This pipeline takes advantage of a genome mapper Bismark, which is capable of alignment DNA sequence data that has been bisulfite treated. Importantly, Bismark wraps around Bowtie2 which does the actual mapping, and provides downstream commands to facilitate removing sequence duplication,and quantify cytosine and thymine coverage at different cytosine motifs across the genome (we focus on CpGs for this study).

## Table of Contents

1. [Data](#data)
2. [Brief Description and Literature on Required Tools and Scripts](#description)
3. [Step 1 - Trimming, adapter removal, and QC](#one)
4. [Step 2 - Creating STAR index](#two)
5. [Step 3 - Mapping with STAR](#three)
6. [Step 4 - Running RSEM](#four)
7. [Step 5 - Filtering, Creating DGEList Object, and Normalization (with limma-voom)](#five)

---

## Data <a name="data"></a>

* [**Link to scripts**](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/tree/master/src/RNAseq)  
* [**Link to data**](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/tree/master/src/MBDBS_seq)
* Reference genome: from NCBI ([GCA_002022765.4 C_virginica-3.0](https://www.ncbi.nlm.nih.gov/genome/?term=crassostrea+virginica))

## Brief Description and Literature on Required Tools and primary R packages <a name="description"></a>

**Trimming and Quality Control**

*TrimGalore!* - a  flexible read trimming tool for Illumina NGS data. Used by dDocent to trim raw RNAseq fragments and remove adapters.

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

Remove reads with poor mapping quality and also cut off 10bp from both the 5' and 3' regions of either strand. This will remove any adapter sequence to improve downstream mapping. This trimming follows the recommendations of bismark when the library prep was down with the pico zymo kit.

* Full Script : []()

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

## Step 2 - Create Reference for Bowtie <a name = "two"></a>

Only need to do this once for all samples (takes a couple of minutes). 

**Command Line**
```
bismark_genome_preparation --bowtie2 --genomic_composition --parallel 10 --verbose /path/toGenomeFolder > bismark_genomePrepartion_log.txt
```

---

Step
