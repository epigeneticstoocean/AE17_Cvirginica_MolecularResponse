
# RNAseq Data QC, Mapping, and Quantification Pipeline

### Overview

This pipeline takes advantage of a genome mapper STAR, which performs transcript alignment by mapping to a reference genome. Importantly, STAR is suited for the de novo discovery of splice junctions, which can be leveraged for identifying novel exons and isoforms. This pipeline couples the STAR mapper with RSEM for transcript quantification. This approach attempts to probabilistically estimate transcript abundance rather than simply count the reads.This may be beneficial for improving transcript count estimates, by probabilistically resolving reads which map to multiple genes (multimappers).  

## Table of Contents

1. [Brief Description and Literature on Required Tools and Scripts](#one)
2. [Step 1 - Creating STAR index](#two)
3. [Step 2 - Mapping with STAR](#three)
4. [Step 3 - Running RSEM](#four)

## Brief Description and Literature on Required Tools and Scripts <a name="one"></a>

**Mapping**

*STAR* - fast RNA-seq aligner than can make to a reference genome and identify identify canonical as well as novel splice junctions. It will output mapped reads as `.sam` or `.bam` files, and with the `--quantMode` it can also create a tab delimited read count output (similar to HT-Seq). In addition, mapped reads can be ouputed as a `.bam` file with transcript coordinates. This can be used downstream by the transcript quantification program RSEM. 

* [Github](https://github.com/alexdobin/STAR)  
* [Publication](https://academic.oup.com/bioinformatics/article/29/1/15/272537)

**Transcript Quantification**

*RSEM* - Transcript quantifier, that can estimate counts at either the transcript (isoform) or gene level. It has a direct workflow with `STAR`, which enables a single line command for both mapping and quantification. Alternatively, it can take`STAR` outputs (specifically `.bam` files with transcript coordinates), and then perform the estimation.

* [Github](https://deweylab.github.io/RSEM/)
* [Publication](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-323)

[**Link to scripts**](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/tree/master/src/RNAseq)
[**Link to data**](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/tree/master/data/RNAseq)

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
