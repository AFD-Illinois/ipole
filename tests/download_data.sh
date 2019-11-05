#!/bin/bash

# Unfortunately you need a BH account to run this, and thus the tests.
# TODO wider distribution...

mkdir test-resources
echo "Copying files from bh.astro.illinois.edu via SSH"
scp bh.astro.illinois.edu:/bd1/ipole-test-resources/*.h5 test-resources/
