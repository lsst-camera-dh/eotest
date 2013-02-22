#!/bin/bash

export PYTHONPATH=.:${PYTHONPATH}

export OUTPUTDIR=../000-00/ccd_read_noise
export CCD_VENDOR=e2v
export SENSOR_ID=000-00

export DB_CREDENTIALS=/nfs/farm/g/lsst/u1/testData/SIMData/pipeline/db_test_app.par

python read_noise_task.py
