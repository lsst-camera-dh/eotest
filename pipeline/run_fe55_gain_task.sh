#!/bin/bash

export PYTHONPATH=.:${PYTHONPATH}

export SENSOR_DIR=/nfs/farm/g/lsst/u1/testData/eotestData/000_00

python fe55_gain_task.py --Vendor e2v --sensor_id 000-00 \
   -d /nfs/farm/g/lsst/u1/testData/SIMData/pipeline/db_test_app.par \
   -f ${SENSOR_DIR}/xray/data/000_00_fe55_0600s_\?\?\[01\].fits -v
