#!/bin/bash

# TODO download if necessary?

for folder in */
do
  folder=${folder: : -1}
  if [[ "$folder" != "test-resources" ]]; then
    cd $folder
    if [[ $folder == *"sample_dump"* ]]; then
      make -f ../../makefile -j8 MODEL=iharm
    else
      make -f ../../makefile -j8 MODEL=thin_disk
    fi
    ./ipole -par $folder.par &> out_$folder.txt
    h5diff -d 1.e-10 image.h5 ../test-resources/$folder.h5
    h5diff -p 1.e-6 image.h5 ../test-resources/$folder.h5
    cd -
  fi
done
