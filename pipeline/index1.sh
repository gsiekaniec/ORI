#!/bin/bash

source ./config.sh

echo "Compute bloom filters..."
howdesbt makebfQ --k=${K} --qgram=${QGRAM_FILE} --bits=${BITS} ${FASTA}
ls *.bf > ./leafname
ORI.py length -g ${DIR} -o ./length.txt