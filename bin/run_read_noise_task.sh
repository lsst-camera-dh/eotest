#!/bin/bash

source ./pipeline_setup.sh

#
# Pipeline example
#
python read_noise_task.py \
    -B ${SENSOR_ID}_BIAS.txt \
    -N ${SENSOR_ID}_SYSNOISE.txt \
    -d ${DB_CREDENTIALS} \
    -s ${SENSOR_ID} -V e2v -o ${SENSOR_ID}/results/read_noise -v

##
## Interactive example
##
#python read_noise_task.py \
#    -b ${DATADIR}/sensorData/${SENSOR_ID}/fe55/debug/${SENSOR_ID}_fe55_bias\*.fits \
#    -n ${DATADIR}/system/noise/debug/noise_\*_debug.fits \
#    -g ${SENSOR_ID}_gains.fits \
#    -s ${SENSOR_ID} -V e2v -o ${SENSOR_ID}/results/read_noise -v
