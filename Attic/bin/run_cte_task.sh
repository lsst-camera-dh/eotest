#!/bin/bash

source ./pipeline_setup.sh

#
# Pipeline example
#
python cte_task.py \
    -F ${SENSOR_ID}_SUPERFLAT.txt \
    -d ${DB_CREDENTIALS} \
    -s ${SENSOR_ID} -V e2v -o ${SENSOR_ID}/results/cte -v
#
##
## Interactive example
##
#python cte_task.py \
#    -f ${DATADIR}/sensorData/${SENSOR_ID}/superflat_500/debug/${SENSOR_ID}_superflat_500_\*.fits \
#    -g ${SENSOR_ID}_gains.fits \
#    -s ${SENSOR_ID} -V e2v -o ${SENSOR_ID}/results/cte -v
