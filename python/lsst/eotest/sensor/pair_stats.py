"""
@brief Compute flat pair statistics.  This implementation is based on
Peter D.'s IDL routine pair_stats.pro.  Use of DM stack enables the
handling of masks when computing the various statistical quantities.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
from __future__ import print_function
from __future__ import absolute_import
import os
import numpy as np
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
from .MaskedCCD import MaskedCCD
import lsst.eotest.image_utils as imutils


class PairStats(object):
    def __init__(self, bias_mean, bias_stddev, flat_mean, flat_var,
                 gain, noise):
        self.bias_mean = bias_mean
        self.bias_stddev = bias_stddev
        self.flat_mean = flat_mean
        self.flat_var = flat_var
        self.gain = gain
        self.noise = noise

    def header(self):
        return """  Exp    Bias Mean   Bias RMS    Flat Mean   Flat Var   Gain e-/DN   Noise e-  
 ------ ----------- ----------- ----------- ----------- ----------- -----------"""

    def summary(self, i=0):
        return ("%5i" % i) + self.__repr__()

    def __repr__(self):
        format = 6*"%12.3F"
        return format % (self.bias_mean, self.bias_stddev,
                         self.flat_mean, self.flat_var,
                         self.gain, self.noise)


def pair_stats(ccd1, ccd2, amp, mask_files=(), binsize=1, bias_frame=None):
    if ccd1.md.get('EXPTIME') != ccd2.md.get('EXPTIME'):
        raise RuntimeError("Exposure times for files %s, %s do not match"
                           % (ccd1.imfile, ccd2.imfile))
    #
    # Mean and variance calculations that account for masks (via
    # ccd1.stat_ctrl, which is the same for both MaskedImages).
    #

    def mean(im): return afwMath.makeStatistics(im, afwMath.MEAN,
                                                ccd1.stat_ctrl).getValue()

    def var(im): return afwMath.makeStatistics(im, afwMath.VARIANCE,
                                               ccd1.stat_ctrl).getValue()
    #
    # Extract imaging region for segments of both CCDs.
    #
    image1 = ccd1.unbiased_and_trimmed_image(amp)
    image2 = ccd2.unbiased_and_trimmed_image(amp)
    #
    # Use serial overscan for bias region.
    #
    b1 = ccd1[amp].Factory(ccd1[amp], ccd1.amp_geom.serial_overscan)
    b2 = ccd2[amp].Factory(ccd2[amp], ccd2.amp_geom.serial_overscan)
    if ccd1.imfile == ccd2.imfile:
        # Don't have pairs of flats, so estimate noise and gain
        # from a single frame, ignoring FPN.
        bmean = mean(b1)
        bvar = var(b1)
        fmean = mean(image1)
        fvar = var(image1)
    else:
        bmean = (mean(b1) + mean(b2))/2.
        #
        # Make a deep copy since otherwise the pixel values in image1
        # would be altered in the ratio calculation.
        #
        fratio_im = afwImage.MaskedImageF(image1, True)
        fratio_im /= image2
        fratio = mean(fratio_im)
        image2 *= fratio
        fmean = (mean(image1) + mean(image2))/2.

        fdiff = afwImage.MaskedImageF(image1, True)
        fdiff -= image2
        fvar = var(fdiff)/2.

        bdiff = afwImage.MaskedImageF(b1, True)
        bdiff -= b2
        bvar = var(bdiff)/2.
    # Compute final stats
    gain = fmean/(fvar - bvar)
    bias_rms = np.sqrt(bvar)
    noise = gain*bias_rms
    return PairStats(bmean, bias_rms, fmean, fvar, gain, noise)


if __name__ == '__main__':
    from lsst.eotest.sensor.sim_tools import simulateFlat

    datadir = '/nfs/slac/g/ki/ki18/jchiang/LSST/SensorTests/test_scripts/work/sensorData/000-00/flat/debug'

    def datapath(x): return os.path.join(datadir, x)
    ccd1 = MaskedCCD(datapath('000-00_flat_005.09s_flat1_debug.fits'))
    ccd2 = MaskedCCD(datapath('000-00_flat_005.09s_flat2_debug.fits'))

    for i, amp in enumerate(ccd1.keys()):
        my_pair_stats = pair_stats(ccd1, ccd2, amp)
        if i == 0:
            print(my_pair_stats.header())
        print(my_pair_stats.summary())
