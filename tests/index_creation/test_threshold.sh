#!/bin/bash

python3 ../../ORI.py threshold -m hamming_matrix.tsv -t 0.0002

diff ./threshold_test.png ./threshold.png
