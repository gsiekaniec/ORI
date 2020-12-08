#!/bin/bash

# kmer size
K="15"
# bloom filters size
BITS="0.2G"
# path to seed file
SEED_FILE="./seed.txt"
# ext, fasta/fastq gz or not
FASTA="*.fasta"

echo "Compute bloom filters..."
howdesbt makebfQ --k=${K} --qgram=${SEED_FILE} --bits=${BITS} ${FASTA}
ls *.bf > ./leafname
ORI.py length -g "." -o ./length.txt
echo "Done."
