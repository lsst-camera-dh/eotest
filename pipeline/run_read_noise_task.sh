#!/bin/bash

export PYTHONPATH=.:${PYTHONPATH}

export NUMBIASFILES=2
export BIAS_00=../000-00/Fe55/Fe55_bias_000-00_00.fits
export BIAS_01=../000-00/Fe55/Fe55_bias_000-00_01.fits

export NUMSYSNOISEFILES=2
export SYSNOISE_00=../000-00/system_noise/system_noise_000-00_00.fits
export SYSNOISE_01=../000-00/system_noise/system_noise_000-00_01.fits

export OUTPUTDIR=../000-00/ccd_read_noise
export CCD_VENDOR=e2v
export SENSOR_ID=000-00

export DB_CREDENTIALS=/nfs/farm/g/lsst/u1/testData/SIMData/pipeline/db_test_app.par

python read_noise_task.py
