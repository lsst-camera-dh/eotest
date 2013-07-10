#!/bin/bash

#
# Pipeline example
#
#python trap_task.py \
#    -f /u/gl/jchiang/ki18/LSST/SensorTests/test_scripts/work/sensorData/000-00/trap/debug/000-00_trap_ppump_debug.fits \
#    -d /nfs/farm/g/lsst/u1/testData/SIMData/pipeline/db_test_app.par \
#    -s 000-00 -V e2v -v

#
# Interactive example
#
python trap_task.py \
    -f /u/gl/jchiang/ki18/LSST/SensorTests/test_scripts/work/sensorData/000-00/trap/debug/000-00_trap_ppump_debug.fits \
    -g 000-00_gains.fits \
    -s 000-00 -V e2v -o traps/data -v
