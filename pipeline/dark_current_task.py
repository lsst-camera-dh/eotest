"""
@brief Dark current task: compute 95th percentile dark current in
units of e-/sec/pixel.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import numpy as np
import pyfits
import lsst.afw.image as afwImage
import image_utils as imutils
from MaskedCCD import MaskedCCD
from TaskParser import TaskParser

parser = TaskParser('Compute 95th percentile dark current.')
parser.add_argument('-f', '--dark_files', type=str,
                    help='file pattern for darks')
parser.add_argument('-F', '--dark_file_list', type=str,
                    help='file contain list of dark files')
parser.add_argument('-t', '--temp_tol', default=1.5, type=float,
                    help='temperature tolerance for CCDTEMP among dark files')
args = parser.parse_args()

dark_files = args.files(args.dark_files, args.dark_file_list)
if args.verbose:
    print 'processing files:', dark_files

sensor_id = args.sensor_id
sensor = args.sensor()
gains = args.system_gains()
mask_files = args.mask_files()

imutils.check_temperatures(dark_files, args.temp_tol)

median_images = {}
md = afwImage.readMetadata(dark_files[0], 1)
for amp in imutils.allAmps:
    median_images[amp] = imutils.fits_median(dark_files, imutils.dm_hdu(amp))
medfile = os.path.join(args.output_dir, '%s_dark_current_map.fits' % sensor_id)
imutils.writeFits(median_images, medfile, md)

ccd = MaskedCCD(medfile, mask_files=mask_files)

dark95s = {}
exptime = md.get('EXPTIME')
print "Segment    95 percentile    median"
for amp in imutils.allAmps:
    image = imutils.unbias_and_trim(ccd[amp].getImage())
    mask = imutils.trim(ccd[amp].getMask())
    imarr = image.getArray()
    mskarr = mask.getArray()
    pixels = imarr.reshape(1, imarr.shape[0]*imarr.shape[1])[0]
    masked = mskarr.reshape(1, mskarr.shape[0]*mskarr.shape[1])[0]
    unmasked = [pixels[i] for i in range(len(pixels)) if masked[i] == 0]
    unmasked.sort()
    unmasked = np.array(unmasked)*gains[amp]/exptime
    dark95s[amp] = unmasked[int(len(unmasked)*0.95)]
    sensor.add_seg_result(amp, 'darkCurrent95', dark95s[amp])
    print "%s         %.2e         %.2e" % (imutils.channelIds[amp],
                                            dark95s[amp],
                                            unmasked[len(unmasked)/2])

dark95mean = np.mean(dark95s.values())
print "CCD: mean 95 percentile value =", dark95mean
sensor.add_ccd_result('darkCurrent95mean', dark95mean)
#
# Update header of dark current median image file with dark files used
# and dark95 values.
#
output = pyfits.open(medfile)
for i, dark in enumerate(dark_files):
    output[0].header.update('DARK%02i' % i, os.path.basename(dark))
for amp in imutils.allAmps:
    output[0].header.update('DARK95%s' % imutils.channelIds[amp], dark95s[amp])
output.writeto(medfile, clobber=True)
