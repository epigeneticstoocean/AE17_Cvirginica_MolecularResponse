#!/bin/bash

files= ls /shared_lab/20180226_RNAseq_2017OAExp/DNAm/rawfiles/2017_exp/17*

echo "Files being trimmed: $files"


for i in $( ls /shared_lab/20180226_RNAseq_2017OAExp/DNAm/rawfiles/2017_exp/17*R1*); do
	for j in $( ls /shared_lab/20180226_RNAseq_2017OAExp/DNAm/rawfiles/2017_exp/17*R2*); do
		R1=$( echo $i | rev | cut -d'/' -f 1 | rev | cut -d'_' -f 1)
		R2=$( echo $j | rev | cut -d'/' -f 1 | rev | cut -d'_' -f 1)
		
		if [ "$R1" == "$R2" ]
		then
			echo "Starting with sample $R1"

			trim_galore \
			--paired \
			--clip_r1 10 \
			--clip_r2 10 \
			--three_prime_clip_R1 10 \
			--three_prime_clip_R2 10 \
			--output_dir /shared_lab/20180226_RNAseq_2017OAExp/DNAm/20190719_fastqc_trim_10bp_Cvirginica_MBD \
			--fastqc_args "--outdir /shared_lab/20180226_RNAseq_2017OAExp/DNAm/20190719_trim_galore_files --threads 18" \
			$i \
			$j \
			2> /shared_lab/20180226_RNAseq_2017OAExp/DNAm/20190719_fastqc_trim_10bp_Cvirginica_MBD/stderr.log
		fi
	done
done
	
