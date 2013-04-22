"""
@brief Bright pixels task: Find pixels and columns in a median image
constructed from an ensemble of darks.  The brightness threshold is
specified in nsig of the noise above the mean.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import numpy as np
import pyfits
import lsst.afw.image as afwImage
import image_utils as imutils
from TaskParser import TaskParser
from BrightPixels import BrightPixels

def _writeFits(images, outfile, exptime):
    output = pyfits.HDUList()
    output.append(pyfits.PrimaryHDU())
    output[0].header['EXPTIME'] = exptime
    for amp in imutils.allAmps:
        output.append(pyfits.ImageHDU(data=images[amp].getArray()))
        output[amp].name = 'AMP%s' % imutils.channelIds[amp]
        output[amp].header.update('DETSIZE', imutils.detsize)
        output[amp].header.update('DETSEC', imutils.detsec(amp))
    output.writeto(outfile, clobber=True)

parser = TaskParser('Find bright pixels and columns')
parser.add_argument('-f', '--files', type=str,
                    help='file pattern for darks')
parser.add_argument('-F', '--file_list', type=str,
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

if args.files is None and args.file_list is None:
    parser.parse_args('--help'.split())

dark_files = parser.files(args.files, args.file_list)
    
if args.verbose:
    print "processing files: ", dark_files
sensor_id = args.sensor_id
sensor = parser.sensor()
gains = parser.system_gains()
mask_files = parser.mask_files()

#
# Check tempertures
#
ccd_temps = [afwImage.readMetadata(x, 1).get('CCDTEMP') for x in dark_files]
temp_avg = imutils.mean(ccd_temps)
tol = args.temp_tol
if max(ccd_temps) - temp_avg > tol or temp_avg - min(ccd_temps) > tol:
    raise RuntimeError("Temperature deviations > %s " % tol +
                       "deg C relative to average.")

median_images = {}
exptime = afwImage.readMetadata(dark_files[0], 1).get('EXPTIME')
for amp in imutils.allAmps:
    median_images[amp] = imutils.fits_median(dark_files, imutils.dm_hdu(amp))
medfile = os.path.join(args.output_dir, '%s_median_dark_bp.fits' % sensor_id)
_writeFits(median_images, medfile, exptime)

bright_pixels = BrightPixels(medfile, mask_files=mask_files,
                             ethresh=args.ethresh, colthresh=args.colthresh,
                             mask_plane=args.mask_plane)

outfile = os.path.join(args.output_dir, '%s_bright_pixel_map.fits' % sensor_id)
total_bright_pixels = 0
print "Segment     # bright pixels"
for amp in imutils.allAmps:
    pixels, columns = bright_pixels.generate_mask(amp, gains[amp], outfile)
    count = len(pixels)
    total_bright_pixels += count
    sensor.add_seg_result(amp, 'numBrightPixels', count)
    print "%s          %i" % (imutils.channelIds[amp], count)

print "Total bright pixels:", total_bright_pixels
sensor.add_ccd_result('numBrightPixels', total_bright_pixels)
