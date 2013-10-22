#!/bin/bash

source ./pipeline_setup.sh

#
# Pipeline example
#
python crosstalk_task.py \
    -F ${SENSOR_ID}_SPOT.txt \
    -d ${DB_CREDENTIALS} \
    -s ${SENSOR_ID} -V e2v -o ${SENSOR_ID}/results/xtalk -v

##
## Interactive example
##
#python crosstalk_task.py \
#    -f ${DATADIR}/sensorData/${SENSOR_ID}/spot/debug/${SENSOR_ID}_spot_\?\?_debug.fits \
#    -s ${SENSOR_ID} -V e2v -o ${SENSOR_ID}/results/xtalk -v
