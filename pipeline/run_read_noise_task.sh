#!/bin/bash

export PYTHONPATH=.:${PYTHONPATH}

#
# Pipeline example
#
python read_noise_task.py \
    -B /nfs/farm/g/lsst/u1/testData/SIMData/pipeline/000-00_BIAS.txt \
    -N /nfs/farm/g/lsst/u1/testData/SIMData/pipeline/000-00_SYSNOISE.txt \
    -d /nfs/farm/g/lsst/u1/testData/SIMData/pipeline/db_test_app.par \
    -s 000-00 -V e2v -o read_noise/data -v

#
# Interactive example
#
#python read_noise_task.py \
#    -b /nfs/farm/g/lsst/u1/testData/SIMData/000-00/Fe55/Fe55_bias_000-00_\*.fits \
#    -n /nfs/farm/g/lsst/u1/testData/SIMData/000-00/system_noise/system_noise_000-00_\*.fits \
#    -g 000-00_gain.fits \
#    -s 000-00 -V e2v -o read_noise/data -v
