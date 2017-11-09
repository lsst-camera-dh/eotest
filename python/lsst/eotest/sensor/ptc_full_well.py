"""
@brief Compute full well.  Fit linear and quadratic functions to data,
with constraint at inferred variance at 1000 e-.  Full well obtains
where the maximum fractional deviation between the linear and
quadratic curves is greater than some maximum value, which is 10% as
specified in LCA-128.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import numpy as np
import numpy.linalg as linalg
from image_utils import allAmps, channelIds

try:
    import pylab_plotter as plot
except ImportError:
    plot = None


def linear_fit_fixedpt(x, y, x0, y0):
    dx = x - x0
    dy = y - y0
    slope = sum(dx*dy)/sum(dx*dx)
    return slope


def quadratic_fit_fixedpt(x, y, x0, y0):
    A = np.matrix([[sum((x**2 - x0**2)**2), sum((x-x0)*(x**2 - x0**2))],
                   [sum((x-x0)*(x**2 - x0**2)), sum((x-x0)**2)]])
    B = [sum((y - y0)*(x**2 - x0**2)), sum((y - y0)*(x - x0))]
    return linalg.solve(A, B)


def gain_est(meanDN, varDN, nmin=10, nmax=100):
    results = np.polyfit(meanDN[nmin:nmax], varDN[nmin:nmax], 1)
    gain = 1./results[0]
    return gain


def full_well(ptcfile, amp, gain=None, fracdevmax=0.10,
              make_plot=False, outfile_prefix=None):
    data = np.recfromtxt(ptcfile)
    data = data.transpose()

    exptime = data[0]
    meanDN = data[(amp-1)*2 + 1]
    varDN = data[(amp-1)*2 + 2]

    if gain is None:
        gain = gain_est(meanDN, varDN)
    #
    # Convert from DN to e-.
    #
    meanNe = gain*meanDN
    varNe = gain*varDN
    #
    # Find the reference e- signal level.
    #
    xref = 1000.
    indx = 0
    while meanNe[indx] < xref:
        indx += 1
    #
    # Linear fit local to this point.
    #
    imin = max(0, indx - 3)
    imax = imin + 5
    f = np.poly1d(np.polyfit(meanNe[imin:imax], varNe[imin:imax], 1))
    fref = f(xref)
    #
    # Loop over signal level points greater than this, compute linear
    # and quadratic fit anchored at xref(=1000e-), and compute maximum
    # deviation.
    #
    dvarmax = 0
    imax = indx + 10
    while dvarmax < fracdevmax and imax < len(meanNe)-1:
        xx = meanNe[indx:imax]
        ff = varNe[indx:imax]

        slope = linear_fit_fixedpt(xx, ff, xref, fref)

        def f1(x): return slope*(x - xref) + fref

        quadfit = quadratic_fit_fixedpt(xx, ff, xref, fref)

        def f2(x): return quadfit[0]*(x**2-xref**2) + quadfit[1]*(x-xref) + fref
        #
        # Here the fractional deviation is computed at the current
        # end-point for the data that are fit.  May need to check the
        # extremum of f1-f2. If "deviation *below* the linear curve"
        # is the criterion, then the current implementation is ok.
        #
        #fracdev = lambda x : np.abs(f1(x) - f2(x))/f1(x)

        def fracdev(x): return (f1(x) - f2(x))/f1(x)
        full_well_est = meanNe[imax]
        dvarmax = fracdev(full_well_est)
        imax += 1

    if make_plot and plot is not None:
        win0 = plot.xyplot(meanNe, varNe, xname='mean(e-)', yname='var(e-)')
        plot.xyplot(xx, ff, oplot=1, color='r')
        x = np.linspace(xref, meanNe[imax], 100)
        plot.curve(x, f1(x), oplot=1, lineStyle=':')
        plot.curve(x, f2(x), oplot=1, lineStyle='--')
        plot.vline(full_well_est)
        win0.set_title('Amp %i, full well = %i e-' % (amp, full_well_est))
        if outfile_prefix is not None:
            plot.save(outfile_prefix + '_ptc.png')

        win1 = plot.curve(x, fracdev(x), xname='mean(e-)',
                          yname='fractional deviation from linear fit')
        plot.hline(fracdevmax)
        plot.vline(full_well_est)
        win1.set_title('Amp %i, full well = %i e-' % (amp, full_well_est))
        if outfile_prefix is not None:
            plot.save(outfile_prefix + '_fracdev.png')

    return full_well_est


if __name__ == '__main__':
    ptcfile = 'ptc_results.txt'
    amp = 1
    print full_well(ptcfile, amp, make_plot=True,
                    outfile_prefix='amp%02i' % amp)
    for amp in allAmps():
        try:
            result = '%2i  %i' % (amp, full_well(ptcfile, amp))
            print result
        except:
            pass
