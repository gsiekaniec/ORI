#!/bin/bash

#no effect, only to draw the threshold on the graph 
THRESHOLD="0.0002" 

echo "Compute bf distance..."
howdesbt distance --list=leafname
ORI.py threshold -m hamming_matrix.tsv -t ${THRESHOLD}
echo "Done: threshold.png"
