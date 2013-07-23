"""
@brief Bright pixels task: Find pixels and columns in a median image
constructed from an ensemble of darks.  The bright pixel threshold is
specified via the --ethresh option and is in units of -e per pixel per
second.  The threshold for the number of bright pixels that define a
bright column is specified via the --colthresh option.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import numpy as np
import pyfits
import lsst.afw.image as afwImage
import image_utils as imutils
from TaskParser import TaskParser
from MaskedCCD import MaskedCCD
from BrightPixels import BrightPixels

parser = TaskParser('Find bright pixels and columns')
parser.add_argument('-f', '--dark_files', type=str,
                    help='file pattern for darks')
parser.add_argument('-F', '--dark_file_list', type=str,
                    help='file containing list of dark files')
parser.add_argument('-e', '--ethresh', default=5, type=int,
                    help='bright pixel threshold in e- per pixel per second')
parser.add_argument('-c', '--colthresh', default=20, type=int,
                    help='bright column threshold in # of bright pixels')
parser.add_argument('-p', '--mask_plane', default='BAD', type=str,
                    help='mask plane to be used for output mask file')
parser.add_argument('-t', '--temp_tol', default=1.5, type=float,
                    help='temperature tolerance for CCDTEMP among dark files')
args = parser.parse_args()

dark_files = args.files(args.dark_files, args.dark_file_list)

if args.verbose:
    print "processing files: ", dark_files

sensor_id = args.sensor_id
sensor = args.sensor()
gains = args.system_gains()
mask_files = args.mask_files()

imutils.check_temperatures(dark_files, args.temp_tol)

median_images = {}
for amp in imutils.allAmps:
    median_images[amp] = imutils.fits_median(dark_files,
                                             imutils.dm_hdu(amp))
    medfile = os.path.join(args.output_dir,
                           '%s_median_dark_bp.fits' % sensor_id)
imutils.writeFits(median_images, medfile, dark_files[0])

ccd = MaskedCCD(medfile, mask_files=mask_files)
md = afwImage.readMetadata(dark_files[0], 1)
exptime = ccd.md.get('EXPTIME')
outfile = os.path.join(args.output_dir,
                       '%s_bright_pixel_map.fits' % sensor_id)
total_bright_pixels = 0
print "Segment     # bright pixels"
for amp in imutils.allAmps:
    bright_pixels = BrightPixels(ccd, amp, exptime, gains[amp])
    pixels, columns = bright_pixels.find()
    bright_pixels.generate_mask(outfile)
    count = len(pixels)
    total_bright_pixels += count
    sensor.add_seg_result(amp, 'numBrightPixels', count)
    print "%s          %i" % (imutils.channelIds[amp], count)

print "Total bright pixels:", total_bright_pixels
sensor.add_ccd_result('numBrightPixels', total_bright_pixels)
