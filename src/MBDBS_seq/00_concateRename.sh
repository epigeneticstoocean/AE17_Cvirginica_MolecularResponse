#!/bin/bash 

filename='/shared_lab/20180226_RNAseq_2017OAExp/DNAm/rawfiles/2017_exp/DNA_seqKey_rv.csv'

while read line; do
	# reading each line
	echo "Matching DNAm ID with sample : $line"

	IFS='    ' # four spaces set as delimiter
	read -ra ADDR <<< "$line" # str is read into an array as tokens separated by IFS
	i="${ADDR[0]}" # Save first column
	j="${ADDR[1]}" # Save second column

	echo "Concatenating sample lists:"
	echo ls /shared_lab/20180226_RNAseq_2017OAExp/DNAm/rawfiles/2017_exp/${i}_*R1.fastq*
	echo "Output file:"
	echo "/shared_lab/20180226_RNAseq_2017OAExp/DNAm/rawfiles/2017_exp/${j}_DNAm_R1.fastq.gz"

	cat /shared_lab/20180226_RNAseq_2017OAExp/DNAm/rawfiles/2017_exp/${i}_*R1.fastq.gz > /shared_lab/20180226_RNAseq_2017OAExp/DNAm/rawfiles/2017_exp/${j}_DNAm_R1.fastq.gz
	cat /shared_lab/20180226_RNAseq_2017OAExp/DNAm/rawfiles/2017_exp/${i}_*R2.fastq.gz > /shared_lab/20180226_RNAseq_2017OAExp/DNAm/rawfiles/2017_exp/${j}_DNAm_R2.fastq.gz
	echo ""
	echo ""
done < $filename
