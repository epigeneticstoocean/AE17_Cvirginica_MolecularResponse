!#/bin/bash

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

echo "Select genome file (.fna format, should include entire path)"
read genome

echo "Select gene annotation file (.gtf, should includ entire path)"
read gene_annotation


STAR --runThreadN 32 \
--runMode genomeGenerate \
--genomeDir $base$output \
--genomeFastaFiles $genome \
--sjdbGTFfile $gene_annotation
