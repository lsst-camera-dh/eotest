import os
import numpy as np
import scipy.stats
import scipy.optimize
import pylab
import lsst.afw.detection as afwDetect
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import image_utils as imutils
from fe55_yield import Fe55Yield

def fe55_lines(x, *args):
    """
    Two Gaussian model of Mn K-alpha and K-beta lines for Fe55 tests.
    The ratio of peak locations is fixed at the line energy ratio, and
    the line widths are assumed to be the same.
    """
    k1, m1, s1, k2 = args
    m2 = 6.49/5.889*m1
    s2 = s1
    value = k1*scipy.stats.norm.pdf(x, loc=m1, scale=s1)
    value += k2*scipy.stats.norm.pdf(x, loc=m2, scale=s2)
    return value

class Xrays(object):
    def __init__(self, image_file, amp):
        self.image_file = image_file
        self.amp = amp
        self.ccdtemp = afwImage.readMetadata(image_file, 1).get('CCDTEMP')
        self.fe55_yield = Fe55Yield(self.ccdtemp)
        hdu = imutils.dm_hdu(amp)
        raw_image = afwImage.ImageF(image_file, hdu)
        self.imarr = raw_image.getArray()
        self.image = imutils.trim(raw_image)
        flags = afwMath.MEANCLIP | afwMath.STDEVCLIP
        stats = afwMath.makeStatistics(self.image, flags)
        self.mean = stats.getValue(afwMath.MEANCLIP)
        self.stdev = stats.getValue(afwMath.STDEVCLIP)
        self.footprint_signal = self._footprint_signal_spans
#        self.footprint_signal = self._footprint_signal_bbox
    def _footprint_signal_spans(self, footprint, buff=None):
        spans = footprint.getSpans()
        total = 0
        for span in spans:
            total += sum(self.imarr[span.getY()][span.getX0():span.getX1()+1])
        return total - footprint.getNpix()*self.mean
    def _footprint_signal_bbox(self, footprint, buff=1):
        bbox = footprint.getBBox()
        xmin = max(imutils.imaging.getMinX(), bbox.getMinX() - buff)
        xmax = min(imutils.imaging.getMaxX(), bbox.getMaxX() + buff)
        ymin = max(imutils.imaging.getMinY(), bbox.getMinY() - buff)
        ymax = min(imutils.imaging.getMaxY(), bbox.getMaxY() + buff)
        subarr = self.imarr[ymin:ymax+1, xmin:xmax+1] 
        signal = sum(subarr.flat) - self.mean*len(subarr.flat)
        return signal
    def pixels(self):
        imarr = self.image.getArray()
        ny, nx = imarr.shape
        return imarr.reshape(1, nx*ny)[0]
    def signals(self, nsig=2, max_npix=9, gain_max=10.):
        sigmin = self.fe55_yield.alpha()/gain_max
        threshold = afwDetect.Threshold(self.mean + nsig*self.stdev)
        fpset = afwDetect.FootprintSet(self.image, threshold)
        signals = np.array([self.footprint_signal(fp) for fp in
                            fpset.getFootprints() if fp.getNpix() < max_npix])
        indx = np.where(signals > sigmin)
        return signals[indx]
    def gain(self, nsig=2, max_npix=9, gain_max=10., make_plot=False):
        signals = self.signals(nsig, max_npix, gain_max)
        flags = afwMath.MEDIAN | afwMath.STDEVCLIP
        stats = afwMath.makeStatistics(signals, flags)
        median = stats.getValue(afwMath.MEDIAN)
        stdev = stats.getValue(afwMath.STDEVCLIP)
        xrange = (median - 10*stdev, median + 10*stdev)
        if make_plot:
            pylab.ion()
            fig = pylab.figure(self.amp)
        else:
            pylab.ioff()
            
        hist = pylab.hist(signals, bins=100, range=xrange,
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
        gain = self.fe55_yield.alpha()/kalpha_peak
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
            pylab.axes().set_title('%s, AMP%s' %
                                   (os.path.basename(self.image_file),
                                    imutils.channelIds[amp]))
        return gain, noise

if __name__ == '__main__':
    data_dir = '/nfs/farm/g/lsst/u1/testData/eotestData/000_00/xray/data' 
    infile = os.path.join(data_dir, '000_00_fe55_0600s_000.fits')
    
    make_plot = True

    print 'AMP   gain    noise'
    for amp in imutils.allAmps[:1]:
        xrays = Xrays(infile, amp)
        gain, noise = xrays.gain(make_plot=make_plot)
        print '%s    %.2f    %.2f' % (imutils.channelIds[amp], gain, noise)
