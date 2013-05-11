#!/bin/bash

export PYTHONPATH=.:${PYTHONPATH}

#
# Pipeline example
#
python flat_pair_task.py \
    -F /nfs/farm/g/lsst/u1/testData/SIMData/pipeline/000-00_FLAT.txt \
    -d /nfs/farm/g/lsst/u1/testData/SIMData/pipeline/db_test_app.par \
    -s 000-00 -V e2v -o flat_pair/data -v

##
## Interactive example
##
#python flat_pair_task.py \
#    -f /nfs/farm/g/lsst/u1/testData/eotestData/000_00/flat/data/\*flat\?.fits \
#    -g 000-00_gain.fits \
#    -s 000-00 -V e2v -o flat_pair/data -v
