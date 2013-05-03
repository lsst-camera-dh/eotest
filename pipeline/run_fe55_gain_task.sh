#!/bin/bash

export PYTHONPATH=.:${PYTHONPATH}

export SENSOR_DIR=/nfs/farm/g/lsst/u1/testData/eotestData/000_00

#
# Pipeline example
#
python fe55_gain_task.py \
    -F /nfs/farm/g/lsst/u1/testData/SIMData/pipeline/000-00_FE55.txt \
    -d /nfs/farm/g/lsst/u1/testData/SIMData/pipeline/db_test_app.par \
    -s 000-00 -V e2v -v
#
# Interactive example
#
#python fe55_gain_task.py \
#    -f ${SENSOR_DIR}/xray/data/000_00_fe55_0600s_\?\?\[01\].fits \
#    -s 000-00 -V e2v -v
