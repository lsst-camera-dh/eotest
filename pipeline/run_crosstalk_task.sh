#!/bin/bash

SENSOR_DIR=/u/gl/jchiang/ki18/LSST/SensorTests/test_scripts/work/sensorData/000-00/spot/debug

##
## Pipeline example
##
#python crosstalk_task.py \
#    -f ${SENSOR_DIR}/000-00_spot_\?\?_debug.fits \
#    -d /nfs/farm/g/lsst/u1/testData/SIMData/pipeline/db_test_app.par \
#    -s 000-00 -V e2v -o xtalk/data -v

#
# Interactive example
#
python crosstalk_task.py \
    -f ${SENSOR_DIR}/000-00_spot_\?\?_debug.fits \
    -s 000-00 -V e2v -o xtalk/data -v
