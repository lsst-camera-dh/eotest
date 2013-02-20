#!/bin/bash

export PYTHONPATH=.:${PYTHONPATH}

export DARKS_LIST=darks.txt
export SENSOR_ID=000-00
export CCD_VENDOR=e2v
export DB_CREDENTIALS=/nfs/farm/g/lsst/u1/testData/SIMData/pipeline/db_test_app.par
export OUTPUTDIR=dark_current/data

python dark_current_task.py
