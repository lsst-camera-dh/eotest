"""
@brief Compute the system gains for each segment in a CCD from Fe55
exposures.  The median gain value from among the input files is
identified as the system gain for that segment and written to the db
tables.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import glob
import argparse
import pyfits

import image_utils as imutils
from ccd250_mask import ccd250_mask
from database.SensorDb import SensorDb, NullDbObject

from xray_gain import hdu_gains

parser = argparse.ArgumentParser(description='Compute gain using Fe55 data')
parser.add_argument('-f', '--files', type=str, 
                    help='file pattern for Fe55 files')
parser.add_argument('-l', '--filelist', type=str,
                    help='file containing list of Fe55 files')
parser.add_argument('-d', '--db_credentials', type=str,
                    help='file containing database credentials')
parser.add_argument('-s', '--sensor_id', type=str,
                    help='sensor ID')
parser.add_argument('-V', '--Vendor', type=str,
                    help='CCD vendor (e.g., e2v, ITL)')
parser.add_argument('-m', '--mask_file', default='ccd250_defects', type=str,
                    help='mask file to use')
parser.add_argument('-o', '--output_dir', default='.', type=str,
                    help='output directory')
parser.add_argument('-v', '--verbose', action='store_true', default=False,
                    help='turn verbosity on')
args = parser.parse_args()

if args.files is None and args.filelist is None:
    parser.parse_args('--help'.split())

#
# Input files. If the filelist option is specified, it takes precedence.
#
if args.filelist is not None:
    Fe55_files = [x.strip() for x in open(args.filelist)]
elif args.files is not None:
    pattern = args.files.replace('\\', '')
    if args.verbose:
        print pattern
    target = os.path.join(pattern)
    Fe55_files = glob.glob(target)
Fe55_files.sort()
    
if args.verbose:
    print "processing files: ", Fe55_files

#
# Database handling
#
if args.db_credentials is not None:
    sensor_id = args.sensor_id
    vendor = args.Vendor
    sensorDb = SensorDb(args.db_credentials)
    sensor = sensorDb.getSensor(vendor, sensor_id, add=True)
else:
    sensor = NullDbObject()
        
#
# Mask file to use.
#
if args.mask_file == 'ccd250_defects':
    mask_files = ('ccd250_defects.fits',)
    ccd250_mask(mask_files[0])
elif mask_file is not None:
    mask_files = (args.mask_file,)
else:
    mask_files = ()

#
# Compute gain distributions for each segment.
#
gain_dists = dict([(amp, []) for amp in imutils.allAmps])
for fe55 in Fe55_files:
    if args.verbose:
        print "processing", fe55
    gains = hdu_gains(fe55, mask_files=mask_files)
    for amp in imutils.allAmps:
        gain_dists[amp].append(gains[amp])
        
seg_gains = dict([(amp, imutils.median(gain_dists[amp]))
                  for amp in imutils.allAmps])

#
# Write output to db table and output file.
#
try:
    os.makedirs(args.output_dir)
except OSError:
    pass
outfile = os.path.join(args.output_dir,
                       "%s_gain.fits" % (sensor_id.replace('-', '_')))
output = pyfits.HDUList()
output.append(pyfits.PrimaryHDU())
    
gain_median = imutils.median(seg_gains.values())
sensor.add_ccd_result('gainMedian', gain_median)
print "Median gain among segments:", gain_median
print "Segment    gain"
for amp in imutils.allAmps:
    sensor.add_seg_result(amp, 'gain', seg_gains[amp])
    print "%s         %.4f" % (imutils.channelIds[amp], seg_gains[amp])
    output[0].header.update("GAIN%s" % imutils.channelIds[amp],
                            seg_gains[amp])
output.writeto(outfile, clobber=True)
