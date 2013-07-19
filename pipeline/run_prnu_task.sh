#!/bin/bash

source ./pipeline_setup.sh

#
# Pipeline example
#
python prnu_task.py \
    -F ${SENSOR_ID}_QE_files.txt \
    -c ${DATADIR}/sensorData/${SENSOR_ID}/superflat_500/debug/${SENSOR_ID}_illumation_correction_debug.fits \
    -d ${DB_CREDENTIALS} \
    -s ${SENSOR_ID} -V e2v -o ${SENSOR_ID}/results/prnu -v

##
## Interactive example
##
#python prnu_task.py \
#    -f ${DATADIR}/sensorData/${SENSOR_ID}/lambda/debug/${SENSOR_ID}_lambda_\*.fits \
#    -g ${SENSOR_ID}_gains.fits \
#    -c ${DATADIR}/sensorData/${SENSOR_ID}/superflat_500/debug/${SENSOR_ID}_illumation_correction_debug.fits \
#    -s ${SENSOR_ID} -V e2v -o prnu/data -v
