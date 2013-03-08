#!/bin/bash
# Example running Fe55 gain interactively

export CCD_VENDOR=e2v
export SENSOR_ID=000-00

python ../pipeline/fe55_gain_task.py --files "/nfs/farm/g/lsst/u1/testData/eotestData/000-00/xray/data/000_00_fe55_0600s_*"
