#!/bin/bash
# Example Running read_noise interactively

export CCD_VENDOR=e2v
export SENSOR_ID=000-00

python ../pipeline/read_noise_task.py --bias "/nfs/farm/g/lsst/u1/testData/eotestData/000-00/xray/data/000_00_fe55_bias_*" --noise "/nfs/farm/g/lsst/u1/testData/eotestData/System/noise/data/sysnoise*.fits" --gain "000_00_gain.fits" -i "000-00" -o "outputDir"
