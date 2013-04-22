#!/bin/bash

export PYTHONPATH=.:${PYTHONPATH}

export OUTPUTDIR=bright_pixels/data

python bright_pixels_task.py \
   -F /nfs/farm/g/lsst/u1/testData/SIMData/pipeline/000-00_DARK.txt \
   -d /nfs/farm/g/lsst/u1/testData/SIMData/pipeline/db_test_app.par \
   -s 000-00 -V e2v -o bright_pixels/data -v
