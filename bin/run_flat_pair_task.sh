#!/bin/bash

source ./pipeline_setup.sh

#
# Pipeline example
#
python flat_pair_task.py \
    -F ${SENSOR_ID}_FLAT.txt \
    -d ${DB_CREDENTIALS} \
    -s ${SENSOR_ID} -V e2v -o ${SENSOR_ID}/results/flat_pair -v

##
## Interactive example
##
#python flat_pair_task.py \
#    -f ${DATADIR}/sensorData/${SENSOR_ID}/flat/debug/${SENSOR_ID}_flat_\*s_flat\*.fits \
#    -g ${SENSOR_ID}_gains.fits \
#    -s ${SENSOR_ID} -V e2v -o ${SENSOR_ID}/results/flat_pair -v
