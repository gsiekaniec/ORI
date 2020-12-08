#!/bin/bash
FASTQ="./test.fastq"
QUAL=9
LENGTH=2000
NB_READS=

#ORI.py suppr_bad_reads -fq ${FASTQ} -q ${QUAL} -l ${LENGTH}

BASEN=$(basename -- $FASTQ)
EXT="${BASEN##*.}"
FASTQ_F="${BASEN%.*}_better_than_${QUAL}.${EXT}"


head -n 16000 ${FASTQ_F} > query_${NB_READS}.fq