#!/bin/bash

export PYTHONPATH=.:${PYTHONPATH}

export SENSOR_ID=000-00
export CCD_VENDOR=e2v
export DB_CREDENTIALS=/nfs/farm/g/lsst/u1/testData/SIMData/pipeline/db_test_app.par
export OUTPUTDIR=bright_pixels/data

python bright_pixels_task.py
