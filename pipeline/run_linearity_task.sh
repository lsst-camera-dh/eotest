#!/bin/bash

export PYTHONPATH=.:${PYTHONPATH}

export FLAT_LIST=linearity_flats.txt
export LINEARITY_OUTFILE=linearity_results.txt
export SENSOR_ID=000-00
export CCD_VENDOR=e2v
export DB_CREDENTIALS=/nfs/farm/g/lsst/u1/testData/SIMData/pipeline/db_test_app.par

python linearity_task.py
