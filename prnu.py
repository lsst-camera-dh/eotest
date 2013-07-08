"""
@brief Pixel response non-uniformity calculation.  This is applied to
the unmasked pixels for an entire CCD.  The computed quantities are
the pixel mean, pixel median, and pixel standard deviation.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import numpy as np
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import image_utils as imutils
from MaskedCCD import MaskedCCD

def extract_unmasked_pixels(masked_image, gain):
    subimage = imutils.unbias_and_trim(masked_image)
    imarr = subimage.getImage().getArray()
    maskarr = subimage.getMask().getArray()
    if imarr.shape != maskarr.shape:
        raise RuntimeError("image and mask subarray shapes differ")
    indx = np.where(maskarr == 0)
    return [x*gain for x in imarr[indx].flat]

def prnu(infile, mask_files, gains, correction_image=None):
    ccd = MaskedCCD(infile, mask_files=mask_files)
    active_pixels = []
    for amp in ccd:
        if correction_image is not None:
            correction = afwImage.ImageF(correction_image, imutils.dm_hdu(amp))
            image = ccd[amp].getImage()
            image /= correction
        active_pixels.extend(extract_unmasked_pixels(ccd[amp], gains[amp]))
    active_pixels = np.array(active_pixels, dtype=np.float)
    flags = afwMath.MEDIAN | afwMath.STDEV
    stats = afwMath.makeStatistics(active_pixels, flags)
    pix_median = stats.getValue(afwMath.MEDIAN)
    pix_stdev = stats.getValue(afwMath.STDEV)
    return pix_stdev, pix_median

if __name__ == '__main__':
    from MaskedCCD import Metadata
    infile = 'work/000-00_lambda_0450.0_20130703151405.fits'
    mask_files = ('work/ccd250_defects.fits', )

    md = Metadata('work/000-00_gains.fits', 1)
    gains = dict([(amp, md.get('GAIN%s' % imutils.channelIds[amp]))
                  for amp in imutils.allAmps])

    pix_stdev, pix_median = prnu(infile, mask_files, gains)

    excess_variance = pix_stdev**2 - pix_median
    print "Excess pixel variance/pixel mean:", excess_variance/pix_median
    if excess_variance > 0:
        print "Fractional excess noise: %.2f" \
            % (np.sqrt(excess_variance)/pix_median*100, )
