#!/bin/bash

python3 ../../ORI.py clean_merge -n leafname -r . -o list_number_file.txt

diff ./list_number_file_test.txt ./list_number_file.txt
