#!/bin/bash

../../HowDeSBT_strains/howdesbt distance --list=leafname --threshold=0.0002 --matrix=hamming_matrix.bin --merge

diff ./0_1_test.bf ./0_1_.bf
