#!/bin/bash

source ./pipeline_setup.sh

#
# Pipeline example
#
python psf_task.py \
    -F ${SENSOR_ID}_FE55.txt \
    -d ${DB_CREDENTIALS} \
    -s ${SENSOR_ID} -V e2v -o ${SENSOR_ID}/results/psf -v

##
## Interactive example
##
#python psf_task.py \
#    -f ${DATADIR}/sensorData/${SENSOR_ID}/fe55/debug/${SENSOR_ID}_fe55_fe55_\*.fits \
#    -g ${SENSOR_ID}_gains.fits \
#    -s ${SENSOR_ID} -V e2v -o ${SENSOR_ID}/results/psf -v

