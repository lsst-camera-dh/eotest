#!/bin/bash

source ./pipeline_setup.sh
#
# The following will be installation-dependent (obviously).
SCRIPTDIR=/u/gl/jchiang/ki18/LSST/SensorTests/test_scripts

#
# Pipeline example
#
python qe_task.py \
    -F ${SENSOR_ID}_QE_files.txt \
    -d ${DB_CREDENTIALS} \
    --ccd_cal_file ${SCRIPTDIR}/qe/OD142.csv \
    --int_sph_cal_file ${SCRIPTDIR}/qe/OD143.csv \
    --wavelength_scan_file ${SCRIPTDIR}/qe/WLscan.txt \
    -s ${SENSOR_ID} -V e2v -o ${SENSOR_ID}/results/qe -v

##
## Interactive example
##
#python qe_task.py \
#    -f ${DATADIR}/sensorData/${SENSOR_ID}/lambda/debug/${SENSOR_ID}_lambda_\*_debug.fits \
#    -g 000-00_gains.fits \
#    --ccd_cal_file ${SCRIPTDIR}/qe/OD142.csv \
#    --int_sph_cal_file ${SCRIPTDIR}/qe/OD143.csv \
#    --wavelength_scan_file ${SCRIPTDIR}/qe/WLscan.txt \
#    -s ${SENSOR_ID} -V e2v -o qe/data -v
