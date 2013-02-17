#!/bin/bash

export PYTHONPATH=.:${PYTHONPATH}

export NUMFE55FILES=2
export FE55_00=../000-00/Fe55/Fe55_exp_000-00_00.fits
export FE55_01=../000-00/Fe55/Fe55_exp_000-00_01.fits

export CCD_VENDOR=e2v
export SENSOR_ID=000-00
export DB_CREDENTIALS=/nfs/farm/g/lsst/u1/testData/SIMData/pipeline/db_test_app.par

python fe55_gain_task.py
