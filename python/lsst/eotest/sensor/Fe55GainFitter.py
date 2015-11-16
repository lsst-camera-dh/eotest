#!/usr/bin/env python
"""
@brief Function to fit the Fe55 signal distribution to obtain the
system gain.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import numpy as np
import scipy.stats
import scipy.optimize
import pyfits
import pylab
import pylab_plotter as plot
import lsst.afw.math as afwMath
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

class Fe55GainFitter(object):
    def __init__(self, signals, ccdtemp=-95):
        self.signals = signals
        self.ccdtemp = ccdtemp
        self._compute_stats()
    def fit(self, xrange=None, bins=100, hist_nsig=10, dADU=50):
        if xrange is None:
            self._set_hist_range(dADU, bins, hist_nsig)
        else:
            self.xrange = xrange
        hist = np.histogram(self.signals, bins=bins, range=self.xrange)
        x = (hist[1][1:] + hist[1][:-1])/2.
        y = hist[0]
        ntot = sum(y)
        #
        # Starting values for two Gaussian fit. The relative
        # normalizations are initially set at the expected line ratio
        # of K-alpha/K-beta = 0.88/0.12.  The relative peak locations
        # and relative widths are fixed in fe55_lines(...) above.
        #
        p0 = (ntot*0.88, self.median, self.stdev/2., ntot*0.12)
        self.pars, pcov = scipy.optimize.curve_fit(fe55_lines, x, y, p0=p0)
        
        kalpha_peak, kalpha_sigma = self.pars[1], self.pars[2]
        kalpha_peak_error = np.sqrt(pcov[1][1])
        fe55_yield = Fe55Yield(self.ccdtemp)
        Ne, Ne_error = fe55_yield.alpha()
        self.gain = Ne/kalpha_peak
#        self.gain_error = float(self.gain*np.sqrt((Ne_error/Ne)**2 +
#                                                  (kalpha_peak_error/kalpha_peak)**2))
        self.gain_error = float(self.gain*kalpha_peak_error/kalpha_peak)
            
        return kalpha_peak, kalpha_sigma
    def _compute_stats(self):
        try:
            flags = afwMath.MEDIAN | afwMath.STDEVCLIP
            stats = afwMath.makeStatistics(self.signals.tolist(), flags)
            self.median = stats.getValue(afwMath.MEDIAN)
            self.stdev = stats.getValue(afwMath.STDEVCLIP)
        except:
            self.median = 0
            self.stdev = 0
    def _set_hist_range(self, dADU, bins, hist_nsig):
        """Set the histogram range for the fit to bracket the estimated
        Kalpha and Kbeta peak locations by +/-dADU"""
        # Set range of histogram to include both Kalpha and Kbeta
        # peaks.
        xmin = max(self.median - hist_nsig*self.stdev, 200)
        xmax = min(self.median*1785./1620. + hist_nsig*self.stdev, 2000)
        xrange = xmin, xmax
        # Determine distribution mode and take that as the initial
        # estimate of the location of the Kalpha peak.
        hist = np.histogram(self.signals, bins=bins, range=xrange)
        xpeak = hist[1][np.where(hist[0] == max(hist[0]))][0]
        # Set final xrange.
        self.xrange = max(0, xpeak - dADU), xpeak*1785./1620. + dADU
    def plot(self, xrange=None, interactive=False, bins=100,
             win=None, subplot=(1, 1, 1), figsize=None, add_labels=False,
             frameLabels=False, amp=1, title=''):
        pylab_interactive_state = pylab.isinteractive()
        pylab.interactive(interactive)
        if win is None:
            if frameLabels:
                xlabel = 'Bias Corrected Event Signal (DN)'
                ylabel = 'Entries / bin'
            else:
                xlabel, ylabel = None, None
            win = plot.Window(subplot=subplot, figsize=figsize,
                              xlabel=xlabel, ylabel=ylabel, size='large')
        else:
            win.select_subplot(*subplot)
        if frameLabels:
            bbox = win.axes[-1].get_position()
            points = bbox.get_points()
            points[0] += 0.025
            points[1] += 0.025
            bbox.set_points(points)
            win.axes[-1].set_position(bbox)
        if xrange is not None:
            self.xrange = xrange
        logscale = True
        if max(self.signals) <= 0:
            logscale = False
        try:
            hist = pylab.hist(self.signals, bins=bins, range=self.xrange,
                              histtype='bar', color='b', log=logscale)
            yrange = 1, max(hist[0])*1.5
            plot.setAxis(self.xrange, yrange)
        except:
            return win
        if add_labels:
            pylab.xlabel('Bias Corrected Event Signal (DN)')
            pylab.ylabel('Entries / bin')
        x = (hist[1][1:] + hist[1][:-1])/2.
        xx = np.linspace(x[0], x[-1], 1000)
        pylab.plot(xx, fe55_lines(xx, *self.pars), 'r--', markersize=3,
                   linewidth=1)
        pylab.annotate(("Amp %i\nGain=%.2f e-/DN") % (amp, self.gain),
                       (0.475, 0.8), xycoords='axes fraction',
                       size='x-small')
        pylab.interactive(pylab_interactive_state)
        return win

if __name__ == '__main__':
    import pyfits
    infile = '/u/gl/jchiang/ki18/LSST/SensorTests/eotest/0.0.0.7/work/results/000-00_psf_results_nsig4.fits'
    chiprob_min = 0.1
    results = pyfits.open(infile)
    for hdu in range(1, 17):
        chiprob = results[hdu].data.field('CHIPROB')
        index = np.where(chiprob > chiprob_min)
        dn = np.array(results[hdu].data.field('DN'), dtype=np.float)[index]

        foo = Fe55GainFitter(dn)
        foo.fit()
        if hdu == 1:
            win = foo.plot(interactive=True, subplot=(4, 4, hdu),
                           figsize=(11, 8.5))
        else:
            foo.plot(interactive=True, subplot=(4, 4, hdu), win=win)
 
            
