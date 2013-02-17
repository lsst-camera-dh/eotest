#!/bin/bash

export PYTHONPATH=.:${PYTHONPATH}

export PTC_FLAT_LIST=ptc_flats.txt
export PTC_OUTFILE=ptc_results.txt
export CCD_VENDOR=e2v
export SENSOR_ID=000-00
export DB_CREDENTIALS=/nfs/farm/g/lsst/u1/testData/SIMData/pipeline/db_test_app.par

python ptc_task.py
