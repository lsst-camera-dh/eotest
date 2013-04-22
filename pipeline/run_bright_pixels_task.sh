#!/bin/bash

export PYTHONPATH=.:${PYTHONPATH}

#
# Pipeline example
#
python bright_pixels_task.py \
   -F /nfs/farm/g/lsst/u1/testData/SIMData/pipeline/000-00_DARK.txt \
   -d /nfs/farm/g/lsst/u1/testData/SIMData/pipeline/db_test_app.par \
   -s 000-00 -V e2v -o bright_pixels/data -v

#
# Interactive example
#
#python bright_pixels_task.py \
#    -f /nfs/farm/g/lsst/u1/testData/eotestData/000_00/dark/data/dark100_\?\?\?.fits \
#    -g 000-00_gain.fits \
#    -s 000-00 -V e2v -o bright_pixels/data -v
