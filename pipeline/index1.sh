#!/bin/bash

# kmer size
K="15"
# bloom filters size
BITS="0.5G"
# path to directory that contains fasta files
DIR="./tests"
# path to seed file
SEED_FILE="${DIR}/seed.txt"
# ext, fasta/fastq gz or not
FASTA="*.fasta"


echo "Compute bloom filters..."
howdesbt makebfQ --k=${K} --qgram=${SEED_FILE} --bits=${BITS} ${DIR}/${FASTA}
ls *.bf > ./leafname
ORI.py length -g ${DIR} -o ./length.txt
