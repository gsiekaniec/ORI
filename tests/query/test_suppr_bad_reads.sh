#!/bin/bash

../../ORI.py suppr_bad_reads -fq ./reads.fastq -q 9 -l 2000

diff ./reads_threshold_test.fastq ./reads_better_than_9.fastq
