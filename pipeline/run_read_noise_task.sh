#!/bin/bash

export PYTHONPATH=.:${PYTHONPATH}

python read_noise_task.py \
    -B /nfs/farm/g/lsst/u1/testData/SIMData/pipeline/000-00_BIAS.txt \
    -N /nfs/farm/g/lsst/u1/testData/SIMData/pipeline/000-00_SYSNOISE.txt \
    -d /nfs/farm/g/lsst/u1/testData/SIMData/pipeline/db_test_app.par \
    -s 000-00 -V e2v -o read_noise/data -v
