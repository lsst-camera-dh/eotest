#!/bin/bash

export PYTHONPATH=.:${PYTHONPATH}

export SENSOR_ID=000-00
export CCD_VENDOR=e2v
export DB_CREDENTIALS=/nfs/farm/g/lsst/u1/testData/SIMData/pipeline/db_test_app.par
export OUTPUTDIR=dark_current/data

#
# Pipeline example
#
python dark_current_task.py \
    -F /nfs/farm/g/lsst/u1/testData/SIMData/pipeline/000-00_DARK.txt \
    -d /nfs/farm/g/lsst/u1/testData/SIMData/pipeline/db_test_app.par \
    -s 000-00 -V e2v -o dark_curr/data -v
#
# Interactive example
#
#python dark_current_task.py \
#    -f /nfs/farm/g/lsst/u1/testData/eotestData/000_00/dark/data/dark100_\?\?\?.fits \
#    -g 000-00_gain.fits \
#    -s 000-00 -V e2v -o dark_curr/data -v


