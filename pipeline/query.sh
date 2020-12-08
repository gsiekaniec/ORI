#!/bin/bash
FASTQ="./JIM8232_1000_reads.fastq"
SEED="./seed.txt"
QUAL=9
LENGTH=2000
NB_READS=4000

echo "Deletion of poor quality reads"
ORI.py suppr_bad_reads -fq ${FASTQ} -q ${QUAL} -l ${LENGTH}

BASEN=$(basename -- $FASTQ)
EXT="${BASEN##*.}"
FASTQ_F="${BASEN%.*}_better_than_${QUAL}.${EXT}"
NUMBER=`expr $NB_READS \* 4`

echo "Get reads..."
head -n $NUMBER ${FASTQ_F} > query_${NB_READS}.fq

echo "Querying index..."
howdesbt queryQ --sort --qgram=$SEED --tree=howde.sbt --threshold=0.5 query_${NB_READS}.fq > results_howde.txt

echo "Matrix creation..."
ORI.py matrix -f results_howde.txt -l leafname -o matrix.tsv

echo "Identification..."
ORI.py identification -m matrix.tsv -f results_howde.txt -le length.txt -l leafname -c $(which clingo) -o ori_results.txt

echo "Formatting results..."
ORI.py beautiful_results -f ori_results.txt -n bf_correspondence.txt --pie_chart	
