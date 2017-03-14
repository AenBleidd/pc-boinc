#!/bin/bash


# clean file of old runs
rm -f boinc_lockfile stderr.txt output/output2.txt

# input data
#   * input/experiments.csv
#   * input/tile.txt

# results stored in
#   * output/
echo bin/pc input/tile2.txt output/output2.txt 0.05 1 2470 0
time bin/pc input/tile2.txt output/output2.txt 0.05 1 2470 0
diff -qs output/output2.txt output/ref_output2.txt
