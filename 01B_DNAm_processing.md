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
6. [Step 6 - Summary Counts by Feature](#six)
6. [Step 7 - DNA Methylation summarized by feature](#seven)

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

* [Github](https://github.com/FelixKrueger/Bismark)

*Bowtie2* - The actual mapper that is used by `Bismark`.

* [Github](https://github.com/BenLangmead/bowtie2)

**General Utility**

* samtools* - A program for basic `bed` file manipulation.  

* [Github](https://github.com/samtools/samtools)

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

### Output

* You will create a series of cytosine reports that will be stored in the designated output folder. These will be used in the next step to make the count matrices and perform the differential methylation.

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

* [`/02_bismarkMapping.sh`](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/src/MBDBS_seq/02_bismarkMapping.sh) : Script that performing mapping,depublication,sorting,and methylation calls for all samples
* [`/03_cytosineReport.sh`](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/src/MBDBS_seq/03_cytosineReport.sh) : Script for creating full cytosine reports

**Disclaimer** Scripts here were hard coded and paths will need to be changed to match personal directories.

## Step 4 - Raw Matrices <a name = "four"></a>

### Overview  

This step uses a custom r script and the methylKit packaged to create a matrix for methylated and unmethylated cytosines, as well as a metaData sheet for storing information about each locus.

### Input

* Mapped sample outputs from `Step 3`, specifically those generated by `coverage2cytosine`.
* [`/AE17_RNAmetaData.RData`](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/data/meta/AE17_RNAmetaData.RData) : Sequence metadata .

### Code

* [Script](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/src/MBDBS_seq/04_methylationMatrix_diffMethylation_methylKit.R)

Reading in Cytosine coverage reports
```
myobj <- methRead(as.list(file.list),
                 sample.id=as.list(ID),
                 assembly="t2",
                 treatment=trt, # Preliminary calling of the treatment
                 mincov = 0, # minimum coverage = 0 (filtering happens separately)
                 pipeline='bismarkCytosineReport') # I use outputs directly from bismark
```

Filtering Samples
```
filtered.myobj <- filterByCoverage(myobj,lo.count=5,lo.perc=NULL,
                                hi.count=200,hi.perc=NULL)
```

Normalizing Counts
```
norm.filtered.myobj <- normalizeCoverage(filtered.myobj,method="median") 
```

Uniting sample counts into a single matrix and destranding
```
meth_filt <- unite(norm.filtered.myobj, destrand=TRUE)
```

Differentrial Expression (generic example)
* Chi-squared testing approach : `test=Chisq`
* Minimum difference in methylation among treatments = 50% : `difference=50`
* Critical q.value : `qvalue=0.01`

```
y<-calculateDiffMeth(x,
                       covariates=cov,test="Chisq",mc.cores=15)
  saveRDS(y,paste0(type,"_",comparison,app,".RData"))
  fwrite(y,paste0(type,"_",comparison,app,".csv"),sep = ",")
  # Multicomparison correction and critical % Diff. Threshold
  y_getDiff <- getMethylDiff(y,difference=50,qvalue=0.01)
```


## Step 5 - CpG intersection and and gene feature counts <a name = "five"></a>


### Overview
In this section we first intersect CpGs identified  in `Step 4 ` with various genomic features (feature tracks generated [previously for Venkataraman et al. 2020](https://github.com/epigeneticstoocean/paper-gonad-meth/tree/master/genome-feature-tracks)) using the the script `05A_CpGIntersectionByFeature.sh`. Next, we Summarize methylation and coverage using the script `05B_CpGCountByFeature.sh`.

### Input
* [Genome Track Files](https://github.com/epigeneticstoocean/paper-gonad-meth/tree/master/genome-feature-tracks)
* Files from `step 4`

### Output
* [`/20200130_CpGbyGeneSummary`](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/tree/master/data/MBDBS_seq/20200130_CpGbyGeneSummary) : Folder with `txt` outputs from `/05B_CpGIntersectionByFeature.sh` script. These included files broken down by feature or DML which contain counts of CpGs covered in individual features (i.e. within a single gene) and summary file summarize methylation within single features by treatment and time.

### Code
  
* [`/05A_CpGCountByFeature.sh`](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/src/MBDBS_seq/05A_CpGCountByFeature.sh)
* [`/05B_CpGIntersectionByFeature.sh`](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/src/MBDBS_seq/05B_CpGIntersectionByFeature.sh)

## Step 6 - Summary Counts by Feature <a name = "six"></a>


### Overview
Combines files generated in `Step 5` into a single `RData` file called [`/gene_CpGcoverageSummary.RData`](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/data/MBDBS_seq/20200130_CpGbyGeneSummary/gene_CpGcoverageSummary.RData), which is used to generated Figure 2A.

### Input
* Files from [`/20200130_CpGbyGeneSummary`](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/tree/master/data/MBDBS_seq/20200130_CpGbyGeneSummary) folder.

### Output
* [`/gene_CpGcoverageSummary.RData`](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/data/MBDBS_seq/20200130_CpGbyGeneSummary/gene_CpGcoverageSummary.RData)

### Code
* [`/CpGSummaryByGeneTable.R`](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/src/Accessory/CpGSummaryByGeneTable.R)

## Step 7 - DNA Methylation summarized by feature <a name = "seven"></a>


### Overview 
Script takes the count (coverage) information for all cpgs with at least 5x coverage starting with the methylKit object, subsets these counts by feature and summarizes (takes the median) methylation for each individual by feature and generates a median methylation by feature summary table with associated meta data for each individuals (i.e. treatment and time).

### Input
* [`/AE17_RNAmetaData.RData`](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/data/meta/AE17_RNAmetaData.RData) : Sequence meta data file.
* [`MBDBS_seq/CpGCoordinates.csv`](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/data/MBDBS_seq/CpGCoordinates.csv) : CpG coordinate list for all CpGs with coverage.
* [`/countMatrix_cov5Filtered_medianNormalized_methylCCounts.csv`](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/data/MBDBS_seq/countMatrix_cov5Filtered_medianNormalized_methylCCounts.csv) : Matrix of unmethylated cytosine counts.
* [`/countMatrix_cov5Filtered_medianNormalized_totalCounts.csv`](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/data/MBDBS_seq/countMatrix_cov5Filtered_medianNormalized_totalCounts.csv) : Matrix of total cytosine counts.
* [`/20200130_IndividualSummaries`](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/tree/master/data/MBDBS_seq/20200130_IndividualSummaries) : CpG counts by feature files.

### Output
* Creates `20200202_AllCountsList_cov5_byFeature.RData` as output. This was reduced to [`/20200202_geneBeta_cov5_byFeature.RData`](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/results/DNAm/20200202_geneBeta_cov5_byFeature.RData) (code is located in [`/Analyses/AE17_fig3_DNAm.R`](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/src/Analyses/AE17_fig3_DNAm.R)).

### Code
* [`/MedianMethylationByFeatureSummary.R`](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/src/Accessory/MedianMethylationByFeatureSummary.R)

