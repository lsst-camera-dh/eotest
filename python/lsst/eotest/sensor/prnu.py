"""
@brief Pixel response non-uniformity calculation.  This is applied to
the unmasked pixels for an entire CCD.  The computed quantities are
the pixel mean, pixel mean, and pixel standard deviation.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
from __future__ import print_function
from __future__ import absolute_import
import os
import numpy as np
import lsst.eotest.image_utils as imutils
from .MaskedCCD import MaskedCCD
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath


def extract_unmasked_pixels(ccd, amp, gain, correction_image=None):
    subimage = ccd.unbiased_and_trimmed_image(amp)
    imarr = subimage.getImage().getArray()
    if correction_image is not None:
        correction = afwImage.ImageF(correction_image, imutils.dm_hdu(amp))
        correction = correction.Factory(correction, ccd.amp_geom.imaging)
        imarr *= correction.getArray()
    maskarr = subimage.getMask().getArray()
    if imarr.shape != maskarr.shape:
        raise RuntimeError("image and mask subarray shapes differ")
    indx = np.where(maskarr == 0)
    return [x*gain for x in imarr[indx].flat]


def prnu(infile, mask_files, gains, bias_frame=None, correction_image=None):
    ccd = MaskedCCD(infile, mask_files=mask_files, bias_frame=bias_frame)
    active_pixels = []
    for amp in ccd:
        active_pixels.extend(extract_unmasked_pixels(ccd, amp, gains[amp],
                                                     correction_image))
    active_pixels = np.array(active_pixels, dtype=float)
    flags = afwMath.MEAN | afwMath.STDEV
    stats = afwMath.makeStatistics(active_pixels, flags)
    pix_mean = stats.getValue(afwMath.MEAN)
    pix_stdev = stats.getValue(afwMath.STDEV)
    return pix_stdev, pix_mean


if __name__ == '__main__':
    infile = 'work/sensorData/000-00/lambda/debug/000-00_lambda_0450.0_debug.fits'
    mask_files = ()

    md = imutils.Metadata('work/000-00_gains.fits', 1)
    gains = dict([(amp, md.get('GAIN%s' % imutils.channelIds[amp]))
                  for amp in imutils.allAmps(infile)])

    pix_stdev, pix_mean = prnu(infile, mask_files, gains)

    excess_variance = pix_stdev**2 - pix_mean
    print("Excess pixel variance/pixel mean:", excess_variance/pix_mean)
    if excess_variance > 0:
        print("Fractional excess noise: %.2f" \
            % (np.sqrt(excess_variance)/pix_mean*100, ))
