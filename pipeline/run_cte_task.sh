#!/bin/bash

export PYTHONPATH=.:${PYTHONPATH}

##
## Pipeline example
##
#python cte_task.py \
#    -F /nfs/farm/g/lsst/u1/testData/SIMData/pipeline/000-00_SUPERFLAT.txt \
#    -d /nfs/farm/g/lsst/u1/testData/SIMData/pipeline/db_test_app.par \
#    -s 000-00 -V e2v -o cte/data -v

#
# Interactive example
#
python cte_task.py \
    -f /home/jchiang/work/LSST/Camera/Sensors/test_scripts/work/sensorData/000-00/superflat_500/130705-104715/000-00_superflat_500_\*.fits \
    -g 000-00_gain.fits \
    -s 000-00 -V e2v -o superflat/data -v
