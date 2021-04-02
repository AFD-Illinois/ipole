#!/bin/bash

# Run all the ipole tests
# Tests have the form test_name/test_name.par
# Different compilation options are usually not necessary,
# except thin_disk, special-cased here

# TODO download the reference files if necessary?

for folder in */
do
  folder=$(echo $folder | sed 's/.$//')
  if [[ "$folder" != "test-resources" ]]; then
    cd $folder
    if [[ $folder == *"sample_dump"* ]]; then
      make -f ../../makefile -j8 MODEL=iharm
    else
      make -f ../../makefile -j8 MODEL=thin_disk
    fi
    ./ipole -par $folder.par &> out_$folder.txt
    ../verify.sh
    cd -
  fi
done
