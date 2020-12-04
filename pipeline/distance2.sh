#!/bin/bash

source ./config.sh

howdesbt distance --list=leafname --threshold=${MERGE_THRESHOLD}
#ORI.py threshold_selection TODO