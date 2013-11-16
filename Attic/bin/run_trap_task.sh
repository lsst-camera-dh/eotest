#!/bin/bash

source ./pipeline_setup.sh

#
# Pipeline example
#
python trap_task.py \
    -f /u/gl/jchiang/ki18/LSST/SensorTests/test_scripts/work/sensorData/${SENSOR_ID}/trap/debug/${SENSOR_ID}_trap_ppump_debug.fits \
    -d ${DB_CREDENTIALS} \
    -s ${SENSOR_ID} -V e2v -o ${SENSOR_ID}/results/traps -v
#
##
## Interactive example
##
#python trap_task.py \
#    -f ${DATADIR}/sensorData/${SENSOR_ID}/trap/debug/${SENSOR_ID}_trap_ppump_debug.fits \
#    -g ${SENSOR_ID}_gains.fits \
#    -s ${SENSOR_ID} -V e2v -o ${SENSOR_ID}/results/traps -v
