#!/bin/bash


# clean file of old runs
rm -f boinc_lockfile stderr.txt

# input data
#   * input/experiments.csv
#   * input/tile.txt

# results stored in
#   * output/
echo bin/pc input/tile.txt output/output.txt 0.05 1 393
time bin/pc input/tile.txt output/output.txt 0.05 1 393
diff -qs output/output.txt output/ref_output.txt
