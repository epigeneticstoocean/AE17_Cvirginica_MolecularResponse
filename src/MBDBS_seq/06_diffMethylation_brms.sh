l#!/bin/bash
set -e
set -u
set -o pipefail


## Use these two lines of code if you have a conda environmnent you'd like to activate
## and list what's in that environment
# source activate msprime-env
# conda list

# This line is necessary for changing the lib path where the correct compiliers are located.
export LD_LIBRARY_PATH=/usr/local/lib:/usr/lib:/usr/local/lib64:/usr/lib64

DIR="/shared_lab/20180226_RNAseq_2017OAExp/DNAm" # set your path
cd $DIR

INPUT="processed_samples/05_countSummary"
RSCRIPT="scripts/R/"
OUTPUT="processed_samples/07_brmsSummary"
DIAG="/diagnostic/"

# Code not error proof - the difference between START and END should be greater than then number of cores
START=200001
END=397063
NCORES=60
range=$(($END-$START+1)) # The range is inclusive (START=1 and END=10 includes both loci 1 and 10)
echo $range
lociPerCore=$(($range/$NCORES))
echo $lociPerCore
start=$(seq $START $lociPerCore $END)
echo $start

for i in $start
do
	echo -e "\n\n"
	end=$((i+$lociPerCore-1))

	if [ $end -gt $END ]
	then
		end=$END
	fi

    ### Call R script
    echo "Running R scripts for lines" echo 'start row' $i 'end row' $end
    #echo ${RSCRIPT}07_diffMethylation_brms.R

    Rscript --vanilla ${RSCRIPT}07_diffMethylation_brms.R ${i} $end ${INPUT} ${OUTPUT} > ${OUTPUT}${DIAG}${i}"_R.out" 2> ${OUTPUT}${DIAG}${i}"_R.error" & echo $!

    # the above line, using the $! at the end, will run the script in the background on a single core
    # in the next iteration, it will run the script in the background on another core for a different set of rows
    # do not wait until the last background process is finished

done

echo "Run complete"
