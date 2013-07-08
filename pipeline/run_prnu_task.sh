#
# Interactive example
#
python prnu_task.py \
    -f ../work/sensorData/000-00/lambda/130703-151308/000-00_lambda_\*.fits \
    -g ../work/000-00_gains.fits \
    -c ../work/sensorData/000-00/superflat_500/130705-104715/000-00_illumation_correction_20130708131311.fits \
    -s 000-00 -V e2v -o prnu/data -v
