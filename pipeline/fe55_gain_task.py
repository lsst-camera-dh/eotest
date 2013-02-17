import os
import sys
import glob

import lsst.afw.math as afwMath

from xray_gain import hdu_gains
from pipeline.file_handling import get_file_list
from database.SensorDb import SensorDb, NullDbObject

median = lambda x : afwMath.makeStatistics(x, afwMath.MEDIAN).getValue()

if __name__ == '__main__':
    if len(sys.argv) == 3:
        sensor_id = sys.argv[1]
        pattern = sys.argv[2].replace('\\', '')
        target = os.path.join(pattern)
        Fe55_files = glob.glob(target)
        Fe55_files.sort()
        pipeline_task = False
    else:
        try:
            Fe55_files = get_file_list('FE55')
            sensor_id = os.environ['SENSOR_ID']
            pipeline_task = True
        except:
            print "usage: python fe55_gain_task.py <sensor_id> <glob pattern>"
            sys.exit(1)

    try:
        vendor = os.environ['CCD_VENDOR']
    except KeyError:
        vendor = 'e2v'

    if pipeline_task:
        sensorDb = SensorDb(os.environ["DB_CREDENTIALS"])
        sensor = sensorDb.getSensor(vendor, sensor_id, add=True)
    else:
        sensor = NullDbObject(vendor, sensor_id)

    nhdu = 16
    
    gain_dists = [[] for x in range(nhdu)]
    for fe55 in Fe55_files:
        print "processing", fe55
        gains = hdu_gains(fe55)
        for hdu in range(nhdu):
            gain_dists[hdu].append(gains[hdu])
    seg_gains = [median(gain_dists[hdu]) for hdu in range(nhdu)]
    
    sensor.add_ccd_result('gainMedian', median(seg_gains))
    print "Median gain among segments:", median(seg_gains)
    print "Segment    gain"
    for hdu in range(nhdu):
        sensor.add_seg_result(hdu, 'gain', seg_gains[hdu])
        print "%02o         %.4f" % (hdu, seg_gains[hdu])
