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

def fe55_gain_fitter(signals, ccdtemp=-95, make_plot=False, xrange=None,
                     bins=100, hist_nsig=10, title='', plot_filename=None):
    """
    Function to fit the distribution of charge cluster DN values from
    a Fe55 dataset.  A two Gaussian model of Mn K-alpha and K-beta
    lines is assumed with the ratio between the K-alpha and K-beta
    energies fixed at 5.889/6.49 and the the Gaussian width of the
    lines set equal.

    The gain (Ne/DN), location and sigma of the K-alpha peak (in units
    of DN) are returned as a tuple.

    If make_plot=True, then a matplotlib plot of the distribution and
    fit is displayed.

    If xrange is not None, then that 2-element tuple is used as the
    histogram x-range.

    If xrange is None, then the histogram x-range is set to 
    +/- hist_nsig*clipped_stdev about the median of the signal
    distribution.
    """
    flags = afwMath.MEDIAN | afwMath.STDEVCLIP
    try:
        stats = afwMath.makeStatistics(signals.tolist(), flags)
    except:
        print signals
        raise
    median = stats.getValue(afwMath.MEDIAN)
    stdev = stats.getValue(afwMath.STDEVCLIP)
    if xrange is None:
        # Set range of histogram to include both Kalpha and Kbeta peaks.
        xmin = max(median - hist_nsig*stdev, 200)
        xmax = min(median*1785./1620. + hist_nsig*stdev, 1000)
        xrange = xmin, xmax
    # Save pylab interactive state.
    pylab_interactive_state = pylab.isinteractive()
    # Determine distribution mode and take that as the location of the
    # Kalpha peak
    hist = np.histogram(signals, bins=bins, range=xrange)
    xpeak = hist[1][np.where(hist[0] == max(hist[0]))][0]
    xrange = max(0, xpeak-200), xpeak*1785./1620. + 200
    if make_plot:
        pylab.ion()
        fig = pylab.figure()
        axes = fig.add_subplot(111)
        hist = pylab.hist(signals, bins=bins, range=xrange,
                          histtype='bar', color='b')
    else:
        pylab.ioff()
        hist = np.histogram(signals, bins=bins, range=xrange)
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
        
    kalpha_peak, kalpha_sigma = pars[1], pars[2]
    fe55_yield = Fe55Yield(ccdtemp)
    gain = fe55_yield.alpha()/kalpha_peak

    if make_plot:
        pylab.xlabel('Bias Corrected Event Signal (DN)')
        pylab.ylabel('Entries / bin')
        xx = np.linspace(x[0], x[-1], 1000)
        pylab.plot(xx, fe55_lines(xx, *pars), 'r--', markersize=3,
                   linewidth=1)
        pylab.annotate(("K-alpha peak = %i DN\n\n" + 
                        "Gain = %.2f e-/DN\n\n")
                       % (kalpha_peak, gain),
                       (0.5, 0.7), xycoords='axes fraction')
        axes.set_title(title)
        if plot_filename is not None:
            pylab.savefig(plot_filename)
    # Restore pylab interactive state.
    pylab.interactive(pylab_interactive_state)
    return gain, kalpha_peak, kalpha_sigma

if __name__ == '__main__':
    import os
    import argparse
    parser= argparse.ArgumentParser(description='System gain from Fe55 data')
    parser.add_argument('psf_par_file',
                        help='Fitted parameters from analysis of Fe55 images')
    parser.add_argument('-a', '--amplifier', type=int, default=1,
                        help='Amplifier (1-16)')
    parser.add_argument('-c', '--chiprob_min', type=float, default=0.1,
                        help='Minimum chi-square probability for cluster fit')
    parser.add_argument('-f', '--fp_est', action='store_true', default=False,
                        help=('Use sum over cluster footprint to estimate DN '
                              + 'instead of 2D Gaussian fit'))
    parser.add_argument('-p', '--plot', action='store_true', default=False,
                        help='Plot distribution and fit')
    parser.add_argument('-o', '--outfile', type=str, default='fe55_dist.png',
                        help='Output file name of plot')
   
    args = parser.parse_args()
    
    results = pyfits.open(args.psf_par_file)
    hdu = args.amplifier
    if args.fp_est:
        dn = np.array(results[hdu].data.field('DN_FP_SUM'), dtype=np.float)
    else:
        dn = np.array(results[hdu].data.field('DN'), dtype=np.float)
    chiprob = results[hdu].data.field('CHIPROB')
    
    indx = np.where(chiprob > args.chiprob_min)

    plot_title = '%s, amp %i' % (os.path.basename(args.psf_par_file),
                                 args.amplifier)
    gain, kalpha_peak, kalpha_sigma = fe55_gain_fitter(dn[indx],
                                                       make_plot=args.plot,
                                                       title=plot_title)
    print "gain for amplifier %i: %s" % (args.amplifier, gain)

    if args.plot:
        pylab.savefig(args.outfile)
