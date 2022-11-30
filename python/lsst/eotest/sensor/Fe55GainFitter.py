#!/usr/bin/env python
"""
@brief Function to fit the Fe55 signal distribution to obtain the
system gain.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
from __future__ import absolute_import
import numpy as np
import scipy.stats
import scipy.optimize
import sklearn.mixture
import pylab
from . import pylab_plotter as plot
import lsst.afw.math as afwMath
from .fe55_yield import Fe55Yield


class Fe55Lines:
    def __init__(self, ccdtemp):
        fe55_yield = Fe55Yield(ccdtemp)
        self.alpha_yield = fe55_yield.alpha()
        self.beta_yield = fe55_yield.beta()

    def __call__(self, x, *args):
        """
        Two Gaussian model of Mn K-alpha and K-beta lines for Fe55 tests.
        The ratio of peak locations is fixed at the line energy ratio and
        the line widths are assumed to be the same.  The branching_ratio
        is the ratio of Kalpha to Kbeta x-rays.
        """
        k1, m1, s1, branching_ratio = args
        k2 = k1/branching_ratio
        m2 = 6.49/5.889*m1
        s2 = s1
        value = k1*scipy.stats.norm.pdf(x, loc=m1, scale=s1)
        value += k2*scipy.stats.norm.pdf(x, loc=m2, scale=s2)
        return value


def fit_gmm(x, n_components=2):
    """
    Fit an n_components Gaussian mixture model to sample data x.

    Return an n_components length sequence of tuples containing
    (weight, mean, sigma) for each Gaussian component, order by
    weight, descending.
    """
    if len(x) < n_components:
        # Use values correponding to a gain of unity.
        return ((0.88, 1590, 40), (0.12, 1590*6.49/5.889, 40))
    gmm = sklearn.mixture.GaussianMixture(n_components=n_components,
                                          warm_start=True)
    gmm = gmm.fit(x[:, np.newaxis])
    gmm = gmm.fit(x[:, np.newaxis])
    gmm = gmm.fit(x[:, np.newaxis])
    component_pars = []
    for weight, mean, covar in zip(gmm.weights_.ravel(), gmm.means_.ravel(),
                                   gmm.covariances_.ravel()):
        component_pars.append((weight, mean, np.sqrt(covar)))
    return sorted(component_pars, key=lambda x: x[0], reverse=True)


class Fe55GainFitter(object):
    def __init__(self, signals, ccdtemp=-95):
        self.signals = signals
        self.fe55_lines = Fe55Lines(ccdtemp)
        self._compute_stats()

    def fit(self, xrange=None, bins=50, hist_nsig=10, dADU=150,
            ythresh_frac=0.1):
        if xrange is None:
            self._set_hist_range(dADU, bins, hist_nsig)
        else:
            self.xrange = xrange
        if self.xrange is not None:
            index = np.where((self.xrange[0] < self.signals) &
                             (self.signals < self.xrange[1]))
            signals = self.signals[index]
        else:
            signals = self.signals

        hist = np.histogram(signals, bins=bins, range=self.xrange)
        x = (hist[1][1:] + hist[1][:-1])/2.
        dx = hist[1][1] - hist[1][0]
        y = hist[0]
        ntot = sum(y)
        index = np.where(y > max(y)*ythresh_frac)
        x = x[index]
        y = y[index]
        # Use sklearn.mixture.GaussianMixture.fit to estimate the
        # Gaussian function parameters for the K-alpha line in the
        # Fe55 cluster DN distribution.
        weight, mean, sigma = fit_gmm(signals, n_components=2)[0]

        # Starting values for two Gaussian fit. The K-alpha mean,
        # sigma, and relative normalizations of the K-alpha and K-beta
        # lines are initially set at the values found from the
        # Gaussian mixture model 2 component fit.  The relative peak
        # location, widths, and relative yields are fixed
        p0 = (weight*ntot*dx, mean, sigma, 7.)

        # Put bounds on the K-alpha peak mean and sigma to prevent the
        # fit from wandering off in parameter space.
        bounds = ((0, 0.99*mean, 0.99*sigma, 4),
                  (np.inf, 1.01*mean, 1.2*sigma, 8))
        try:
            self.pars, pcov = scipy.optimize.curve_fit(self.fe55_lines, x, y,
                                                       p0=p0, bounds=bounds)
        except:
            self.pars = p0
            pcov = [[0, 0], [0, 0]]

        kalpha_peak, kalpha_sigma = self.pars[1], self.pars[2]
        kalpha_peak_error = np.sqrt(pcov[1][1])
        Ne, Ne_error = self.fe55_lines.alpha_yield
        self.gain = Ne/kalpha_peak
        self.gain_error \
            = self.gain*np.sqrt((Ne_error/Ne)**2 +
                                (kalpha_peak_error/kalpha_peak)**2)

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
        xmin = self.median - hist_nsig*self.stdev
        xmax = self.median*1785./1620. + hist_nsig*self.stdev
        xrange = xmin, xmax
        # Determine distribution mode and take that as the initial
        # estimate of the location of the Kalpha peak.
        hist = np.histogram(self.signals, bins=bins, range=xrange)
        xpeak = hist[1][np.where(hist[0] == max(hist[0]))][0]
        # Set final xrange.
        self.xrange = max(0, xpeak - dADU), xpeak*1785./1620. + dADU
        if self.xrange[1] <= self.xrange[0]:
            self.xrange = None

    def plot(self, xrange=None, interactive=False, bins=50,
             win=None, subplot=(1, 1, 1), figsize=None, add_labels=False,
             frameLabels=False, amp=1, title='', xrange_scale=1):
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
        # Expand x-axis scale.
        xrange_mid = sum(self.xrange)/2.
        xmin = (min(self.xrange) - xrange_mid)*xrange_scale + xrange_mid
        xmax = (max(self.xrange) - xrange_mid)*xrange_scale + xrange_mid
        logscale = True
        if max(self.signals) <= 0:
            logscale = False
        try:
            hist = pylab.hist(self.signals, bins=bins, range=(xmin, xmax),
                              histtype='bar', color='b', log=logscale)
            yrange = 1, max(hist[0])*1.5
            plot.setAxis((xmin, xmax), yrange)
        except:
            return win
        if add_labels:
            pylab.xlabel('Bias Corrected Event Signal (DN)')
            pylab.ylabel('Entries / bin')
        x = (hist[1][1:] + hist[1][:-1])/2.
        xx = np.linspace(x[0], x[-1], 1000)
        ydata = xrange_scale*self.fe55_lines(xx, *self.pars)
        pylab.plot(xx, ydata, 'r--', markersize=3, linewidth=1)
        pylab.annotate(("Amp %i\nGain=%.2f e-/DN") % (amp, self.gain),
                       (0.475, 0.8), xycoords='axes fraction',
                       size='x-small')
        pylab.interactive(pylab_interactive_state)
        return win


if __name__ == '__main__':
    import astropy.io.fits as fits
    infile = '/u/gl/jchiang/ki18/LSST/SensorTests/eotest/0.0.0.7/work/results/000-00_psf_results_nsig4.fits'
    chiprob_min = 0.1
    results = fits.open(infile)
    for hdu in range(1, 17):
        chiprob = results[hdu].data.field('CHIPROB')
        index = np.where(chiprob > chiprob_min)
        dn = np.array(results[hdu].data.field('DN'), dtype=float)[index]

        foo = Fe55GainFitter(dn)
        foo.fit()
        if hdu == 1:
            win = foo.plot(interactive=True, subplot=(4, 4, hdu),
                           figsize=(11, 8.5))
        else:
            foo.plot(interactive=True, subplot=(4, 4, hdu), win=win)
