#!/bin/bash

MERGE_THRESHOLD="0.0002"

echo "Merge bf..."
howdesbt distance --list=leafname --matrix=hamming_matrix.bin --threshold=${MERGE_THRESHOLD} --merge
ORI.py clean_merge -n path/to/leafname -r "." -o ./bf_correspondence.txt

ls *.bf > leafname
ORI.py merge_length -b ./leafname -l ./length.txt -c ./bf_correspondence.txt -o ./merge_length.txt
mv ./merge_length.txt ./length.txt
echo "Done."
