#!/bin/bash

source ./pipeline_setup.sh

#
# Pipeline example
#
python fe55_gain_task.py \
    -F ${SENSOR_ID}_FE55.txt \
    -d ${DB_CREDENTIALS} \
    -s ${SENSOR_ID} -V e2v -o ${SENSOR_ID}/results/fe55 -v

##
## Interactive example
##
#python fe55_gain_task.py \
#    -f ${DATADIR}/sensorData/${SENSOR_ID}/fe55/debug/${SENSOR_ID}_fe55_fe55\*.fits \
#    -s ${SENSOR_ID} -V e2v -O ${SENSOR_ID}_gains.fits -v
