"""
@brief Useful functions for pipeline tasks.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
from __future__ import print_function
import os
import subprocess
from database.SensorDb import SensorDb, NullDbObject
from database.SensorGains import SensorGains


def setVariable(key, value):
    subprocess.call('pipelineSet %s %s' % (key, value), shell=True)


def get_file_list(prefix, sensor_id):
    infile = '%s_%s.txt' % (sensor_id, prefix)
    print("Reading filenames from %s" % infile)
    my_files = [x.strip() for x in open(infile) if x[0] != '#']
    return my_files


def export_file_list(files, prefix, sensor_id):
    outfile = '%s_%s.txt' % (sensor_id, prefix)
    print("Write filenames to %s" % outfile)
    output = open(outfile, 'w')
    for item in files:
        output.write('%s\n' % item)
    output.close()


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
            print("Setting system gain to 5.5 e-/DN for all segments.")
            gains = SensorGains(5.5)
    return gains, sensor
