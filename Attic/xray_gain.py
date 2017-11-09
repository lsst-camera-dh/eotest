"""
@brief Re-implementation of P. Doherty's IDL function to compute
sensor gain from Fe55 data.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import numpy as np
import lsst.afw.math as afwMath
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDetect
import image_utils as imutils
from MaskedCCD import MaskedCCD
from fe55_yield import Fe55Yield


class Fe55Gain(object):
    def __init__(self, imfile, mask_files=(), ccdtemp_par=-100):
        self.ccd = MaskedCCD(imfile, mask_files=mask_files)
        self.md = afwImage.readMetadata(imfile, 0)
        try:
            ccdtemp = self.md.get('CCDTEMP')
        except:
            ccdtemp = ccdtemp_par
        self.fe55_yield = Fe55Yield(ccdtemp).alpha()
        self._footprint_signal = self._footprint_signal_spans
        #self._footprint_signal = self._footprint_signal_bbox

    def _generate_stats(self, amp):
        #
        # Extract the imaging region of the masked image for the
        # specified amplifier.
        #
        image = imutils.trim(self.ccd[amp])
        #
        # Store numpy array of full segment for
        # self._footprint_signal, since footprint spans use full
        # segment image pixel coordinates.
        #
        self.arr = self.ccd[amp].getImage().getArray()
        stats = afwMath.makeStatistics(image,
                                       afwMath.STDEVCLIP | afwMath.MEDIAN,
                                       self.ccd.stat_ctrl)
        self.noise = stats.getValue(afwMath.STDEVCLIP)
        self.median = stats.getValue(afwMath.MEDIAN)

    def _footprint_signal_spans(self, footprint, buff=None):
        spans = footprint.getSpans()
        total = 0
        for span in spans:
            total += sum(self.arr[span.getY()][span.getX0():span.getX1()+1])
        return total - footprint.getNpix()*self.median

    def _footprint_signal_bbox(self, footprint, buff=1):
        bbox = footprint.getBBox()
        xmin = max(imutils.imaging.getMinX(), bbox.getMinX() - buff)
        xmax = min(imutils.imaging.getMaxX(), bbox.getMaxX() + buff)
        ymin = max(imutils.imaging.getMinY(), bbox.getMinY() - buff)
        ymax = min(imutils.imaging.getMaxY(), bbox.getMaxY() + buff)
        subarr = self.arr[ymin:ymax+1, xmin:xmax+1]
        signal = sum(subarr.flat) - self.median*len(subarr.flat)
        return signal

    def gain(self, amp, jmargin=None, max_npix=9, buff=1):
        if jmargin is None:
            jmargins = range(10)
        else:
            jmargins = (jmargin,)
        self._generate_stats(amp)
        values = []
        fp_sets = []
        for j in jmargins:
            margin = self.noise*(j + 3)
            threshold = afwDetect.Threshold(self.median + margin)
            fp_set = afwDetect.FootprintSet(self.ccd[amp], threshold)
            fp_sets.append(fp_set)
            signals = [self._footprint_signal(fp, buff) for fp in
                       fp_set.getFootprints() if fp.getNpix() < max_npix]
            try:
                stats = afwMath.makeStatistics(signals, afwMath.MEANCLIP)
                mean = stats.getValue()
                values.append(self.fe55_yield/mean)
            except:
                pass
        my_gain = np.median(values)
        return my_gain


def hdu_gains(infile, mask_files=()):
    fe55 = Fe55Gain(infile, mask_files=mask_files)
    gains = {}
    for amp in imutils.allAmps:
        gains[amp] = fe55.gain(amp)
    return gains


if __name__ == '__main__':
    from simulation.sim_tools import CCD
    from ccd250_mask import ccd250_mask

    mask_file = 'CCD250_DEFECTS_mask.fits'
    ccd250_mask(mask_file)

    test_file = 'fe55_test_image.fits'

    ccd = CCD(exptime=1, gain=3, ccdtemp=-100)
    ccd.add_bias(level=2000, sigma=2)
    ccd.add_dark_current(level=2e-3)
    ccd.add_Fe55_hits(nxrays=1000)
    ccd.writeto(test_file)

    gains = hdu_gains(test_file, mask_files=(mask_file,))

    for key, value in gains.items():
        print key, value
