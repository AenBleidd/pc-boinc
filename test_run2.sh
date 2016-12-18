#!/bin/bash


# clean file of old runs
rm -f boinc_lockfile stderr.txt

# input data
#   * input/experiments.csv
#   * input/tile.txt

# results stored in
#   * output/
echo bin/pc input/tile2.txt output/output2.txt 0.05 1 2470
time bin/pc input/tile2.txt output/output2.txt 0.05 1 2470
diff -qs output/output2.txt output/ref_output2.txt
