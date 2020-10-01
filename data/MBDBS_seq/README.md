# MBD-BSseq Data

Summaries of MBD-BSseq data. Files were generated following the pipeline described [`HERE`](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/01B_DNAm_processing.md). Generally, these files were created using `Trim Galore!`, `Bismark`, and `methylKit` for the trimming, aligning, and filtering steps, respectively.

### Folders

* [`/20200130_CpGbyGeneSummary`](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/tree/master/data/MBDBS_seq/20200130_CpGbyGeneSummary) : Feature count summaries
* [`/20200130_IndividualSummaries`](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/tree/master/data/MBDBS_seq/20200130_IndividualSummaries) : DNA methylation of CpGs by feature
* [`bed`](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/tree/master/data/MBDBS_seq/bed) : Bed files.

### Files

* [`/methylKitObj_cov5Filtered_medianNormalized_united.RData`](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/data/MBDBS_seq/methylKitObj_cov5Filtered_medianNormalized_united.RData) : MBD-BSseq summary counts in RData file (all counts).

* [`/countMatrix_cov5Filtered_medianNormalized_methylCCounts.csv`](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/data/MBDBS_seq/countMatrix_cov5Filtered_medianNormalized_methylCCounts.csv) : Methylated Cytosines counts 
* [`/countMatrix_cov5Filtered_medianNormalized_totalCounts.csv`](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/data/MBDBS_seq/countMatrix_cov5Filtered_medianNormalized_totalCounts.csv) : Total Cytosines counts (coverage)
* [`/countMatrix_cov5Filtered_medianNormalized_unMethylCCounts.csv`](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/data/MBDBS_seq/countMatrix_cov5Filtered_medianNormalized_unMethylCCounts.csv) : Unmethylated Cytosines counts
* [`/countMatrix_cov5Filtered_medianNormalized_beta.tar.xz`](https://github.com/epigeneticstoocean/AE17_Cvirginica_MolecularResponse/blob/master/data/MBDBS_seq/countMatrix_cov5Filtered_medianNormalized_beta.tar.xz): DNA Methylation proportion, beta (compressed)
