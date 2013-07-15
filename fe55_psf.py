import numpy as np
from scipy.special import erf, gammaincc
import scipy.optimize
import lsst.afw.detection as afwDetect
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import image_utils as imutils
from MaskedCCD import MaskedCCD
import pylab_plotter as plot

_sqrt2 = np.sqrt(2)

def pixel_integral(x, y, x0, y0, sigma):
    """
    Integrate 2D Gaussian centered at (x0, y0) with width sigma over a
    square pixel at (x, y) with unit width.
    """
    x1, x2 = x - 0.5, x + 0.5
    y1, y2 = y - 0.5, y + 0.5

    Fx = 0.5*(erf((x2 - x0)/_sqrt2/sigma) - erf((x1 - x0)/_sqrt2/sigma))
    Fy = 0.5*(erf((y2 - y0)/_sqrt2/sigma) - erf((y1 - y0)/_sqrt2/sigma))
    
    return Fx*Fy

def psf_func(pos, *args):
    x0, y0, N0, sigma = args
    if type(pos) == type([]):
        return N0*np.array([pixel_integral(x[0], x[1], x0, y0, sigma) 
                            for x in pos])
    return N0*pixel_integral(x[0], x[1], x0, y0, sigma)

def chisq(pos, dn, args):
    return sum((psf_func(pos, *tuple(args)) - np.array(dn))**2)

class PsfGaussFit(object):
    def __init__(self, nsig=3, min_npix=5, gain=5, Ne0=1590, sigma0=0.36):
        self.nsig = nsig
        self.min_npix = min_npix
        self.gain = gain
        self.Ne0 = Ne0
        self.sigma0 = sigma0
        self.sig_fit, self.chiprob = [], []
    def process_image(self, image):
        image -= imutils.bias_image(image)
        try:
            imarr = image.getArray()
        except AttributeError:
            imarr = image.getImage().getArray()

        flags = afwMath.MEDIAN | afwMath.STDEVCLIP
        statistics = afwMath.makeStatistics(image, flags) 
        median = statistics.getValue(afwMath.MEDIAN)
        stdev = statistics.getValue(afwMath.STDEVCLIP)

        threshold = afwDetect.Threshold(median + self.nsig*stdev)
        fpset = afwDetect.FootprintSet(image, threshold)

        for fp in fpset.getFootprints():
            if fp.getNpix() < self.min_npix:
                continue
            spans = fp.getSpans()
            positions = []
            dn = []
            peak = [pk for pk in fp.getPeaks()][0]
            p0 = (pk.getIx(), pk.getIy(), self.Ne0/self.gain, self.sigma0)
            for span in spans:
                y = span.getY()
                for x in range(span.getX0(), span.getX1() + 1):
                    positions.append((x, y))
                    dn.append(imarr[y][x])
            try:
                pars, _ = scipy.optimize.curve_fit(psf_func, positions, 
                                                   dn, p0=p0)
                self.sig_fit.append(pars[3])
                chi2 = chisq(positions, dn, pars)
                dof = fp.getNpix() - 4
                self.chiprob.append(gammaincc(dof/2., chi2/2.))
            except RuntimeError:
                pass
    def results(self, min_prob=0.1):
        sig_fit = np.array(self.sig_fit, dtype=np.float)
        chiprob = np.array(self.chiprob, dtype=np.float)
        indx = np.where(chiprob > min_prob)
        return sig_fit[indx], chiprob[indx]

if __name__ == '__main__':
    #infile = 'simulation/sensorData/000-00/fe55/debug/000-00_fe55_fe55_00_debug.fits'
    infile = 'fe55_0060s_000.fits'

    ccd = MaskedCCD(infile)

    fitter = PsfGaussFit()
    for amp in imutils.allAmps:
        print 'processing amp:', amp
        fitter.process_image(ccd[amp])

    sig_fit, chiprob = fitter.results()

    plot.pylab.ion()

    flags = afwMath.MEDIAN | afwMath.STDEVCLIP
    stats = afwMath.makeStatistics(sig_fit, flags)
    median = stats.getValue(afwMath.MEDIAN)
    stdev = stats.getValue(afwMath.STDEVCLIP)
    plot.histogram(sig_fit, xname='Fitted sigma values',
                   xrange=(median-3*stdev, median+3*stdev))

    plot.xyplot(chiprob, sig_fit, xname='chi-square prob.', yname='sigma')
