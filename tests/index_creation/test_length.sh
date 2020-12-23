#!/bin/bash

python3 ../../ORI.py length -g . -o length.txt -s 15 -fpr 0.05

diff ./bf_min_size_test.txt ./bf_min_size.txt
diff ./length_test.txt ./length.txt
