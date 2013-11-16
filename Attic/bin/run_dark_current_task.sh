#!/bin/bash

source ./pipeline_setup.sh

#
# Pipeline example
#
python dark_current_task.py \
    -F ${SENSOR_ID}_DARK.txt \
    -d ${DB_CREDENTIALS} \
    -s ${SENSOR_ID} -V e2v -o ${SENSOR_ID}/results/dark_curr -v

##
## Interactive example
##
#python dark_current_task.py \
#    -f ${DATADIR}/sensorData/${SENSOR_ID}/dark/debug/${SENSOR_ID}_dark_dark_\*_debug.fits \
#    -g ${SENSOR_ID}_gains.fits \
#    -s ${SENSOR_ID} -V e2v -o ${SENSOR_ID}/results/dark_curr -v
