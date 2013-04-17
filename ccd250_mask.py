import pyfits
import lsst.afw.detection as afwDetect
import lsst.afw.image as afwImage
import lsst.daf.base as dafBase
import image_utils as imutils
from BrightPixels import BrightPixels
from simulation.sim_tools import CCD

mask_image = 'ccd250_mask_image.fits'
mask_file = 'CCD250_DEFECTS_mask.fits'
mask_plane = 'CCD250_DEFECTS'

outer_edge_width = 10
bloom_stop_width = 5

signal = 10

#
# Create an empty set of frames to fill with an image of the mask from
# which we will generate the masks.
#
ccd = CCD(exptime=1)

#
# Write the output file with a primary HDU so that the DMstack code
# can append only image extensions (and not write to the PHDU).
#
hdulist = pyfits.HDUList()
hdulist.append(pyfits.PrimaryHDU())
hdulist.writeto(mask_file, clobber=True)

#
# Amplifiers 1 (AMP10), 8 (AMP17), 9 (AMP07) and 16 (AMP00) are
# along the perimeter.
#
# Amps 8 & 16 have roll-off near the prescan:
# 
xmin = imutils.imaging.getMinX()
xmax = xmin + outer_edge_width
for amp in (8, 16):
    imarr = ccd.segments[amp].image.getArray()
    imarr[:, xmin:xmax] += signal
#
# Amps 1 & 9 have roll-off near the serial overscan:
#
xmax = imutils.imaging.getMaxX() + 1
xmin = xmax - outer_edge_width
for amp in (1, 9):
    imarr = ccd.segments[amp].image.getArray()
    imarr[:, xmin:xmax] += signal
#
# Loop over all amps, set signal in perimeter and around blooming
# stop; write masks to the FITS file.
#
ymax = imutils.imaging.getMaxY() + 1
ymin = ymax - bloom_stop_width
for i, amp in enumerate(imutils.allAmps):
    image = ccd.segments[amp].image
    imarr = image.getArray()
    #
    # Set signal in row direction along perimeter.
    #
    imarr[0:outer_edge_width, :] += signal
    #
    # Set signal around blooming stop
    #
    imarr[ymin:ymax, :] += signal
ccd.writeto(mask_image)
#
# Use BrightPixels code to write the mask file.
#
gain = 1
mask = afwImage.MaskU(image.getDimensions())
mask.addMaskPlane(mask_plane)
for amp in imutils.allAmps:
    bright_pixels = BrightPixels(mask_image, amp)
    pixels, columns = bright_pixels.find(gain, ethresh=signal/2.)
    bright_pixels.write_mask(mask_file, mask_plane)
