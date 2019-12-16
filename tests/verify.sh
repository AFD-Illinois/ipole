#!/bin/bash

# Lists differences between a test run and reference image
# Run this from a specific test's directory

folder=$(basename $PWD)

echo "ABSOLUTE 1e-3"
h5diff -d 1.e-3 image.h5 ../test-resources/$folder.h5 /pol
echo "ABSOLUTE 1e-5"
h5diff -d 1.e-5 image.h5 ../test-resources/$folder.h5 /pol
echo "ABSOLUTE 1e-8"
h5diff -d 1.e-8 image.h5 ../test-resources/$folder.h5 /pol
echo ""
echo ""
echo "RELATIVE 1e-3"
h5diff -p 1.e-3 image.h5 ../test-resources/$folder.h5 /pol
echo "RELATIVE 1e-5"
h5diff -p 1.e-5 image.h5 ../test-resources/$folder.h5 /pol
echo "RELATIVE 1e-8"
h5diff -p 1.e-8 image.h5 ../test-resources/$folder.h5 /pol

python ../../scripts/plot_pol.py image.h5 > plot_img.txt
python ../../scripts/plot_pol.py ../test-resources/$folder.h5 > plot_ref.txt
diff plot_img.txt plot_ref.txt
