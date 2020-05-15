#!/bin/bash

# Run all the ipole tests
# Tests have the form test_name/test_name.par
# Different compilation options are usually not necessary,
# except thin_disk, special-cased here

# TODO download the reference files if necessary?

for folder in */
do
  folder=${folder: : -1}
  if [[ "$folder" != "test-resources" ]]; then
    cd $folder
    if [[ $folder == *"thin_disk"* ]]; then
      make -f ../../makefile -j8 MODEL=thin_disk
    else
      make -f ../../makefile -j8 MODEL=iharm
    fi
    ./ipole -par $folder.par &> out_$folder.txt
    ../verify.sh
    cd -
  fi
done
