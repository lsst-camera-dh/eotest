"""
@brief Useful functions for pipeline tasks.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import subprocess
from database.SensorDb import SensorDb, NullDbObject
from database.SensorGains import SensorGains

def setVariable(key, value):
    subprocess.call('pipelineSet %s %s' % (key, value), shell=True)

def get_file_list(prefix):
    numfiles = int(os.environ["NUM%sFILES" % prefix])
    my_files = []
    for i in range(numfiles):
        my_files.append(os.environ["%s_%02i" % (prefix, i)])
    return my_files

def export_file_list(files, prefix):
    setVariable("NUM%sFILES" % prefix, "%s" % len(files))
    for i, item in enumerate(files):
        setVariable("%s_%02i" % (prefix, i), item)

def setup(argv, indx):
    try:
        sensor_id = os.environ['SENSOR_ID']
        vendor = os.environ['CCD_VENDOR']
        sensorDb = SensorDb(os.environ['DB_CREDENTIALS'])
        sensor = sensorDb.getSensor(vendor, sensor_id)
        gains = SensorGains(vendor=vendor, vendorId=sensor_id)
    except KeyError:
        sensor = NullDbObject()
        try:
            gains = SensorGains(float(argv[indx]))
        except IndexError:
            print "Setting system gain to 5.5 e-/DN for all segments."
            gains = SensorGains(5.5)
    return gains, sensor
