#!/bin/bash

python3 ../../ORI.py merge_length -b leafname_merge -l length.txt -c list_number_file.txt -o merge_length.txt

diff ./merge_length_test.txt ./merge_length.txt
