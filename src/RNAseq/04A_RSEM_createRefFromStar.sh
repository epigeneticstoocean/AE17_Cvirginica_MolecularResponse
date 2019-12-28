!#/bin/bash

rsem-prepare-reference \
--gtf /shared_lab/20180226_RNAseq_2017OAExp/RNA/references/gene_annotation/KM_CV_genome_edit.gtf \
--star \
-p 8 \
/shared_lab/20180226_RNAseq_2017OAExp/RNA/references/genome/GCF_002022765.2_C_virginica-3.0_genomic.fna \
/shared_lab/20180226_RNAseq_2017OAExp/RNA/references/RSEM_ref2/RSEM_ref2

