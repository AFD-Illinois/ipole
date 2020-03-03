#!/bin/bash

# Lists differences between a test run and reference image
# Run this from a specific test's directory

folder=$(basename $PWD)

echo "Pixels different by >1% absolute"
h5diff -d 1.e-2 ../test-resources/${folder}_image.h5 image.h5 /pol
echo "Pixels different by >1% relative"
h5diff -p 1.e-2 ../test-resources/${folder}_image.h5 image.h5 /pol

python ../../scripts/compare.py ../test-resources/${folder}_image.h5 image.h5
