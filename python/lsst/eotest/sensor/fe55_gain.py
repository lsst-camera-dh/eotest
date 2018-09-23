"""
@brief Compute system gain from Fe55 data.  Use afwDetect to find
thresholded footprints of Fe55 hits, and fit the histogram of the ADU
signals associated with those footprints with a pair of Gaussian lines
with mean ratios fixed to the Mn Kalpha/Kbeta line energy ratio.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
from __future__ import print_function
from __future__ import absolute_import
import os
import numpy as np
import scipy.stats
import scipy.optimize
import pylab
import lsst.afw.detection as afwDetect
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.eotest.image_utils as imutils
from .MaskedCCD import MaskedCCD
from .fe55_yield import Fe55Yield


def fe55_lines(x, *args):
    """
    Two Gaussian model of Mn K-alpha and K-beta lines for Fe55 tests.
    The ratio of peak locations is fixed at the line energy ratio, and
    the line widths are assumed to be the same.
    """
    k1, m1, s1, k2 = args
    m2 = 6.49/5.889*m1
#    m2 = 1778./1620.*m1  # These are the nominal Kbeta/Kalpha e- from
    # Janesick, but the Kbeta number looks wrong.
    s2 = s1
    value = k1*scipy.stats.norm.pdf(x, loc=m1, scale=s1)
    value += k2*scipy.stats.norm.pdf(x, loc=m2, scale=s2)
    return value


class Xrays(object):
    _fig_num = 0

    def __init__(self, ccd, amp, bg_reg=(10, 10)):
        self.ccd = ccd
        self.amp = amp
        self.ccdtemp = ccd.md.get('CCDTEMP')
        self.fe55_yield = Fe55Yield(self.ccdtemp)
        raw_image = ccd[amp]
        try:
            self.imarr = raw_image.getArray()
        except AttributeError:
            self.imarr = raw_image.getImage().getArray()
        self.image = imutils.trim(raw_image, imaging=ccd.amp_geom.imaging)
        self.image -= self._bg_image(*bg_reg)
        flags = afwMath.MEANCLIP | afwMath.STDEVCLIP
        stats = afwMath.makeStatistics(self.image, flags, self.ccd.stat_ctrl)
        self.mean = stats.getValue(afwMath.MEANCLIP)
        self.stdev = stats.getValue(afwMath.STDEVCLIP)
        self.footprint_signal = self._footprint_signal_spans
        #self.footprint_signal = self._footprint_signal_bbox

    def _bg_image(self, nx, ny):
        bg_ctrl = afwMath.BackgroundControl(nx, ny, self.ccd.stat_ctrl)
        bg = afwMath.makeBackground(self.ccd[self.amp], bg_ctrl)
        image_region = self.ccd.amp_geom.imaging
        return imutils.trim(bg.getImageF(), imaging=image_region)

    def _footprint_signal_spans(self, footprint, buff=None):
        spans = footprint.getSpans()
        total = 0
        for span in spans:
            total += sum(self.imarr[span.getY()][span.getX0():span.getX1()+1])
        return total - footprint.getArea()*self.mean

    def _footprint_signal_bbox(self, footprint, buff=1):
        imaging = self.ccd.amp_geom.imaging
        bbox = footprint.getBBox()
        xmin = max(imaging.getMinX(), bbox.getMinX() - buff)
        xmax = min(imaging.getMaxX(), bbox.getMaxX() + buff)
        ymin = max(imaging.getMinY(), bbox.getMinY() - buff)
        ymax = min(imaging.getMaxY(), bbox.getMaxY() + buff)
        subarr = self.imarr[ymin:ymax+1, xmin:xmax+1]
        signal = sum(subarr.flat) - self.mean*len(subarr.flat)
        return signal

    def pixels(self):
        imarr = self.image.getArray()
        ny, nx = imarr.shape
        return imarr.reshape(1, nx*ny)[0]

    def signals(self, nsig=2, max_npix=9, gain_max=10.):
        sigmin = self.fe55_yield.alpha()[0]/gain_max
        threshold = afwDetect.Threshold(self.mean + nsig*self.stdev)
        fpset = afwDetect.FootprintSet(self.image, threshold)
        signals = np.array([self.footprint_signal(fp) for fp in
                            fpset.getFootprints() if fp.getArea() < max_npix])
        indx = np.where(signals > sigmin)
        return signals[indx]

    def gain(self, nsig=2, max_npix=9, gain_max=10., make_plot=False,
             xrange=None, bins=100, hist_nsig=10):
        signals = self.signals(nsig, max_npix, gain_max)
        flags = afwMath.MEDIAN | afwMath.STDEVCLIP
        stats = afwMath.makeStatistics(signals, flags)
        median = stats.getValue(afwMath.MEDIAN)
        stdev = stats.getValue(afwMath.STDEVCLIP)
        if xrange is None:
            # Set range of histogram to include both Kalpha and Kbeta peaks.
            xmin = median - hist_nsig*stdev
            xmax = median*1785./1620. + hist_nsig*stdev
            xrange = xmin, xmax
        if make_plot:
            pylab.ion()
            fig = pylab.figure(self._fig_num)
            self._fig_num += 1
        else:
            pylab.ioff()

        hist = pylab.hist(signals, bins=bins, range=xrange,
                          histtype='bar', color='b')
        x = (hist[1][1:] + hist[1][:-1])/2.
        y = hist[0]
        ntot = sum(y)
        #
        # Starting values for two Gaussian fit. The relative
        # normalizations are initially set at the expected line ratio
        # of K-alpha/K-beta = 0.88/0.12.  The relative peak locations
        # and relative widths are fixed in fe55_lines(...) above.
        #
        p0 = (ntot*0.88, median, stdev/2., ntot*0.12)
        pars, _ = scipy.optimize.curve_fit(fe55_lines, x, y, p0=p0)

        kalpha_peak = pars[1]
        gain = self.fe55_yield.alpha()[0]/kalpha_peak
        noise = self.stdev*gain

        if make_plot:
            pylab.xlabel('Bias Corrected Event Signal (DN)')
            pylab.ylabel('Entries / bin')
            xx = np.linspace(x[0], x[-1], 1000)
            pylab.plot(xx, fe55_lines(xx, *pars), 'r--', markersize=3,
                       linewidth=1)
            pylab.annotate(("K-alpha peak = %i DN\n\n" +
                            "Gain = %.2f e-/DN\n\n" +
                            "Noise = %.2f e-")
                           % (kalpha_peak, gain, noise),
                           (0.1, 0.7), xycoords='axes fraction')
        return gain, noise


def hdu_gains(fe55_file, mask_files=()):
    ccd = MaskedCCD(fe55_file, mask_files=mask_files)
    ccd.setAllMasks()
    ccdtemp = ccd.md.get('CCDTEMP')
    gains = {}
    for amp in ccd:
        xrays = Xrays(ccd, amp)
        gains[amp], noise = xrays.gain()
    return gains


if __name__ == '__main__':
    from .MaskedCCD import MaskedCCD

    data_dir = '/nfs/farm/g/lsst/u1/testData/eotestData/000_00/xray/data'
    infile = os.path.join(data_dir, '000_00_fe55_0600s_000.fits')
#    data_dir = '/nfs/slac/g/ki/ki18/jchiang/LSST/SensorTests/test_scripts/work/sensorData/000-00/fe55/debug'
#    infile = os.path.join(data_dir, '000-00_fe55_fe55_00_debug.fits')

    make_plot = True

    ccd = MaskedCCD(infile)

    print('AMP   gain    noise')
    for amp in list(ccd.keys())[:1]:
        xrays = Xrays(ccd, amp)
        gain, noise = xrays.gain(make_plot=make_plot)
        print('%s    %.2f    %.2f' % (imutils.channelIds[amp], gain, noise))
