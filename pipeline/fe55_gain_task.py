"""
@brief Compute the system gains for each segment in a CCD from Fe55
exposures.  The median gain value from among the input files is
identified as the system gain for that segment and written to the db
tables.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import sys
import glob
import argparse
import pyfits

import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath
import image_utils as imUtils
import pipeline.pipeline_utils as pipeUtils
from database.SensorDb import SensorDb, NullDbObject

from xray_gain import hdu_gains

median = lambda x : afwMath.makeStatistics(x, afwMath.MEDIAN).getValue()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(\
                      description='Compute gain using Fe55 data')
    parser.add_argument('-f', '--files', help="input Fe55 file pattern", type=str)
    parser.add_argument('-v', '--verbose', help="turn verbosity on", \
             action='store_true', default=False)
    args = parser.parse_args()

    if args.files:
        pattern = args.files.replace('\\', '')
        target = os.path.join(pattern)
        Fe55_files = glob.glob(target)
        Fe55_files.sort()
        pipeline_task = False
        if args.verbose:
            print "processing files: ", Fe55_files
    else:
        try:
            sensor_id = os.environ['SENSOR_ID']
            Fe55_files = pipeUtils.get_file_list('FE55', sensor_id)
            pipeline_task = True
        except:
            print "usage: python fe55_gain_task.py <Fe55 file pattern>"
            sys.exit(1)

    try:
        sensor_id = os.environ['SENSOR_ID']
        vendor = os.environ['CCD_VENDOR']
        sensorDb = SensorDb(os.environ["DB_CREDENTIALS"])
        sensor = sensorDb.getSensor(vendor, sensor_id, add=True)
    except KeyError:
        sensor = NullDbObject()

    # Read the full segment.
    bbox = afwGeom.Box2I()
    
#    # Omit the first 200 columns to avoid edge roll-off.
#    bbox = afwGeom.Box2I(afwGeom.Point2I(200, 0),
#                         afwGeom.Point2I(541, 2021))

    gain_dists = dict([(amp, []) for amp in imUtils.allAmps])
    for fe55 in Fe55_files:
        print "processing", fe55
        gains = hdu_gains(fe55, bbox=bbox)
        for amp in imUtils.allAmps:
            gain_dists[amp].append(gains[amp])
    seg_gains = [imUtils.median(gain_dists[amp]) for amp in imUtils.allAmps]

    if not pipeline_task:
        outfile = "%s_gain.fits" % (sensor_id.replace('-', '_'))
        output = pyfits.HDUList()
        output.append(pyfits.PrimaryHDU())
    
    
    sensor.add_ccd_result('gainMedian', imUtils.median(seg_gains))
    print "Median gain among segments:", imUtils.median(seg_gains)
    print "Segment    gain"
    for amp in imUtils.allAmps:
        sensor.add_seg_result(amp, 'gain', seg_gains[amp-1])
        print "%s         %.4f" % (imUtils.channelIds[amp], seg_gains[amp-1])
        if not pipeline_task:
            output[0].header.update("GAIN%s" % imUtils.channelIds[amp], seg_gains[amp-1])
  
    if not pipeline_task:
        output.writeto(outfile, clobber=True)
