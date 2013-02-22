#!/bin/bash

export PYTHONPATH=.:${PYTHONPATH}

export CCD_VENDOR=e2v
export SENSOR_ID=000-00
export DB_CREDENTIALS=/nfs/farm/g/lsst/u1/testData/SIMData/pipeline/db_test_app.par

python fe55_gain_task.py
