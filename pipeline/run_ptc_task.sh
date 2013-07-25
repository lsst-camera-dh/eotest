#!/bin/bash

export PYTHONPATH=.:${PYTHONPATH}

export OUTPUTDIR=ptc/data
export CCD_VENDOR=e2v
export SENSOR_ID=000-00
export DB_CREDENTIALS=/nfs/farm/g/lsst/u1/testData/SIMData/pipeline/db_test_app.par

python ptc_task.py
