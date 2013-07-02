#!/bin/bash

export PYTHONPATH=.:${PYTHONPATH}

#
# Pipeline example
#
#python qe_task.py \
#   -F /nfs/farm/g/lsst/u1/testData/SIMData/pipeline/000-00_QE_files.txt \
#   -d /nfs/farm/g/lsst/u1/testData/SIMData/pipeline/db_test_app.par \
#   -s 000-00 -V e2v -o qe/data -v

#
# Interactive example
#
python qe_task.py \
    -f /nfs/farm/g/lsst/u1/testData/HarvardData/112-01/final/bss70/qe/112_01_qe_\[0-9\]\*.fits.gz\
    -g 000-00_gains.fits \
    --ccd_cal_file /u/gl/jchiang/ki18/LSST/SensorTests/test_scripts/qe/OD142.csv \
    --int_sph_cal_file /u/gl/jchiang/ki18/LSST/SensorTests/test_scripts/qe/OD143.csv \
    --wavelength_scan_file /u/gl/jchiang/ki18/LSST/SensorTests/test_scripts/qe/WLscan.txt \
    -s 000-00 -V e2v -o qe/data -v
