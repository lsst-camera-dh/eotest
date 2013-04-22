"""
@brief Compute the system gains for each segment in a CCD from Fe55
exposures.  The median gain value from among the input files is
identified as the system gain for that segment and written to the db
tables.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import glob
import pyfits
import image_utils as imutils
from TaskParser import TaskParser
from xray_gain import hdu_gains

parser = TaskParser('Compute gain using Fe55 data')
parser.add_argument('-f', '--files', type=str, 
                    help='file pattern for Fe55 files')
parser.add_argument('-F', '--file_list', type=str,
                    help='file containing list of Fe55 files')
args = parser.parse_args()

if args.files is None and args.file_list is None:
    parser.parse_args('--help'.split())

#
# Input files. If the filelist option is specified, it takes precedence.
#
Fe55_files = parser.files(args.files, args.file_list)
    
if args.verbose:
    print "processing files: ", Fe55_files

sensor_id = args.sensor_id
sensor = parser.sensor()
mask_files = parser.mask_files()        

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
