"""
@brief Class to compute detector response characteristics (e.g., full well,
linearity) from flat pairs data.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import sys
import numpy as np
import astropy.io.fits as fits
import scipy.optimize
import lsst.eotest.image_utils as imutils
import pylab
from . import pylab_plotter as plot


def _fwc_solve(f1_pars, f2_pars, g=0.1):
    """
    The solution (provided by e2v) for the flux corresponding to the
    full-well capacity described in LCA-10103-A.  This is simply the
    flux at which the quadratic fit to the data in the vicinity of
    full well lies below the linear extrapolation of the detector
    response from lower flux levels by a fraction g.
    """
    a, b, c = f2_pars
    d, f = f1_pars
    x = (-np.sqrt((b + d*g - d)**2 - 4.*a*(c + f*g - f)) - b - d*g + d)/2./a
    return x


class DetectorResponse(object):
    def __init__(self, infile, ptc=None, gain_range=None,
                 hdu_name='DETECTOR_RESPONSE'):
        if infile[-5:] == '.fits':
            self._read_from_fits(infile, hdu_name)
        else:
            self._read_from_text(infile)
        self._sort_by_flux()
        self._index = {}
        #self._compute_gain_selection(ptc, gain_range, infile)

    def _compute_gain_selection(self, ptc, gain_range, infile):
        self._index = {}
        if ptc is None or gain_range is None:
            return
        if len(ptc[1].data.field('AMP01_MEAN')) != len(self.flux):
            raise RuntimeError('Number of measurements in PTC file ' +
                               'differs from detector response file.')
        for amp in imutils.allAmps(infile):
            mean = ptc[1].data.field('AMP%02i_MEAN' % amp)
            var = ptc[1].data.field('AMP%02i_VAR' % amp)
            gain = mean/var
            self._index[amp] = np.where((gain >= gain_range[0]) &
                                        (gain <= gain_range[1]))

    def _sort_by_flux(self):
        index = np.argsort(self.flux)
        self.flux = self.flux[index]
        for amp in self.Ne:
            self.Ne[amp] = self.Ne[amp][index]

    def _read_from_fits(self, infile, hdu_name='DETECTOR_RESPONSE'):
        all_amps = imutils.allAmps(infile)
        with fits.open(infile) as foo:
            hdu = foo[hdu_name]
            self.flux = np.fabs(np.array(hdu.data.field('FLUX'),
                                         dtype=np.float))
            self.Ne = dict([(amp, np.array(hdu.data.field('AMP%02i_SIGNAL' % amp),
                                           dtype=np.float)) for amp in all_amps])

    def _read_from_text(self, infile):
        data = np.recfromtxt(infile)
        data = data.transpose()
        self.flux = data[0]
        self.Ne = dict([(amp, ne) for amp, ne in
                        zip(imutils.allAmps(), data[1:])])

    def full_well(self, amp, max_non_linearity=0.02,
                  frac_offset=0.1, make_plot=False, plotter=None,
                  multipanel=False, fit_range=(1e2, 5e4)):
        if plotter is None:
            plotter = plot
        max_frac_dev, f1_pars, Ne, flux = self.linearity(amp,
                                                         fit_range=fit_range)
        f1 = np.poly1d(f1_pars)
        dNfrac = 1 - Ne/f1(flux)
        indexes = np.arange(len(dNfrac))
        good_vals = np.where(np.abs(dNfrac) <= max_non_linearity)[0]
        if good_vals.sum() < 2:
            print("Not enough good points to fit full_well %i %i"
                  % (amp, good_vals.sum()))
            return (0., f1)

        imin = good_vals[-1]
        try:
            imax = np.where((np.abs(dNfrac) >= frac_offset) &
                            (indexes > imin))[0][0] + 1
        except IndexError:
            imax = len(dNfrac)
        #imin = imax - 3  # this is the proposed change from e2v on 2015-06-04
        x, y = flux[imin:imax], Ne[imin:imax]
        f2_pars = np.polyfit(x, y, 2)
        f2 = np.poly1d(f2_pars)
        fwc = _fwc_solve(f1_pars, f2_pars)
        full_well_est = f2(fwc)
        # Save pylab interactive state.
        pylab_interactive_state = plotter.pylab.isinteractive()
        Ne_scale = 1e-3
        if make_plot:
            if multipanel:
                plotter.xyplot(flux, Ne*Ne_scale, xname='', yname='',
                               new_win=False)
            else:
                plotter.xyplot(flux, Ne,
                               xname='illumination (flux*exptime, arb. units)',
                               yname='Ne (e- per pixel)')
            plotter.xyplot(flux[imin:imax], Ne[imin:imax]*Ne_scale, oplot=1,
                           color='r', markersize=6)
            xx = np.linspace(min(x), max(x), 50)
            plotter.curve(flux, f1(flux)*Ne_scale, oplot=1, color='g')
            plotter.curve(xx, f2(xx)*Ne_scale, oplot=1, color='b')
            plotter.hline(full_well_est*Ne_scale)
            plotter.vline(fwc)
        # Restore pylab interactive state.
        plotter.pylab.interactive(pylab_interactive_state)
        return full_well_est, f1

    def full_well_polyfit(self, amp, order=15, fit_range=(1e2, 5e4),
                          poly_fit_min=1e3, frac_offset=0.1, dfrac=None,
                          make_plot=False, plotter=None, multipanel=False):
        if plotter is None:
            plotter = plot
        flux = self.flux
        Ne = self.Ne[amp]
        if self._index:
            flux = flux[self._index[amp]]
            Ne = Ne[self._index[amp]]
        #
        # Fit linear part of response curve.
        #
        max_Ne_index = np.where(Ne == max(Ne))[0][0]
        indx1 = np.where((Ne > fit_range[0]) & (Ne < fit_range[1])
                         & (flux <= flux[max_Ne_index]))
        f1 = np.poly1d(np.polyfit(flux[indx1], Ne[indx1], 1))
        #
        # Fit polynomial of specified order.
        #
        Ne_vals = f1(flux)
        dNe_frac = (Ne - Ne_vals)/Ne_vals
        if dNe_frac[-1] > -frac_offset:
            message = """Detector response data do not extend past the
              specified fractional offset of %.1f%% for the full
              well definition.""" % (frac_offset*100,)
            raise RuntimeError(message)
        #
        # Solve for specified fractional difference between
        # linear and polynomial fits to determine full well.
        #
        if dfrac is not None:
            indxp = np.where((dNe_frac > -(frac_offset+dfrac)) &
                             (Ne > poly_fit_min))
        else:
            indxp = np.where((Ne > poly_fit_min))
        fp = np.poly1d(np.polyfit(flux[indxp], Ne[indxp], order))

        def df(xx): return 1 - fp(xx)/f1(xx) - frac_offset
        x = flux[indxp]
        imin = np.where(x > 1e1)[0][0]
        flux0 = scipy.optimize.brentq(df, x[imin], x[-1])
        full_well = int(fp(flux0))

        # Save pylab interactive state.
        pylab_interactive_state = plotter.pylab.isinteractive()
        if make_plot:
            yrange = (0, max(max(Ne), full_well)*1.1)
            plotter.xyplot(flux, Ne,
                           xname='Illumination (flux*exptime, arb. units)',
                           yname='e- per pixel',
                           yrange=yrange, new_win=(not multipanel))
            plotter.curve(flux, f1(flux), oplot=1, color='g')
            plotter.xyplot(flux[indxp], Ne[indxp], oplot=1, color='r')
            xx = np.linspace(0, max(flux), 100)
            plotter.curve(xx, fp(xx), oplot=1, color='b')
            plotter.vline(flux0)
            plotter.hline(full_well)
            plotter.pylab.annotate('Amplifier %s\nfull well = %i'
                                   % (amp, int(full_well)),
                                   (0.1, 0.8), xycoords='axes fraction')
        # Restore pylab interactive state.
        plotter.pylab.interactive(pylab_interactive_state)
        return full_well, fp

    def plot_diagnostics(self, flux, Ne, indxp, f1, fp):
        plot.pylab.ion()
        plot.xyplot(flux, Ne)
        plot.xyplot(flux[indxp], Ne[indxp], oplot=1, color='r')
        plot.curve(flux, f1(flux), oplot=1)
        plot.curve(flux, fp(flux), oplot=1, color='b')

    def linearity(self, amp, fit_range=None, spec_range=(1e3, 9e4)):
        flux, Ne = self.flux, self.Ne[amp]
        if self._index:
            # Apply selection to remove points with outlier gains from
            # mean-variance estimate.
            flux = flux[self._index[amp]]
            Ne = Ne[self._index[amp]]
        if fit_range is None:
            fit_range = spec_range
        max_Ne_index = np.where(Ne == max(Ne))[0][0]
        index = np.where((Ne > fit_range[0]) & (Ne < fit_range[1])
                         & (flux <= flux[max_Ne_index]))
        if sum(index[0]) < 1:
            print(f"No selected points to fit linearity for amp {amp}")
            return (0., (1, 0), Ne, flux)
        # Fit a linear slope to these data, using the variance for the
        # signal levels assuming Poisson statistics in the chi-square
        # and fixing the y-intercept to zero.  Package the slope as
        # part of a tuple to be passed to np.poly1d.
        slope = len(Ne[index])/np.sum(flux[index]/Ne[index])
        f1_pars = slope, 0
        f1 = np.poly1d(f1_pars)
        # Further select points that are within the specification range
        # for computing the maximum fractional deviation.
        spec_index = np.where((Ne > spec_range[0]) & (Ne < spec_range[1])
                             & (flux <= flux[max_Ne_index]))

        flux_spec = flux[spec_index]
        Ne_spec = Ne[spec_index]
        if len(Ne_spec) < 1:
            print(f"No selected points for max frac dev for amp {amp}")
            return (0., (1, 0), Ne, flux)

        dNfrac = 1 - Ne_spec/f1(flux_spec)
        return max(abs(dNfrac)), f1_pars, Ne, flux, Ne[index], flux[index]

    def plot_linearity(self, maxdev, f1_pars, Ne, flux, max_dev=0.02):
        top_rect = [0.1, 0.3, 0.8, 0.6]
        bottom_rect = [0.1, 0.1, 0.8, 0.2]
        fig = pylab.figure()
        top_ax = fig.add_axes(top_rect)
        bot_ax = fig.add_axes(bottom_rect, sharex=top_ax)

        f1 = np.poly1d(f1_pars)
        dNfrac = 1 - Ne/f1(flux)
        # Plot flux vs e-/pixel.
        top_ax.loglog(Ne, flux, 'ko')
        top_ax.loglog(Ne, f1(Ne), 'r-')
        top_ax.set_ylabel('flux')
        for label in top_ax.get_xticklabels():
            label.set_visible(False)

        # Plot fractional residuals vs e-/pixel.
        bot_ax.semilogx(Ne, dNfrac, 'ko')
        bot_ax.semilogx(Ne, np.zeros(len(Ne)), 'r-')
        plot.setAxis(yrange=(-1.5*max_dev, 1.5*max_dev))
        bot_ax.set_ylabel('fractional residual flux')
        bot_ax.set_xlabel('e-/pixel')


if __name__ == '__main__':
    infile = '000-00_det_response.txt'
    detResp = DetectorResponse(infile)
    make_plot = False
    print("Amp       max. dev.   full well")
    for amp in imutils.allAmps():
        max_dev, fit_pars = detResp.linearity(amp, make_plot=make_plot)
        sys.stdout.write("%2i         %.3f       " % (amp, max_dev))
        try:
            full_well = detResp.full_well(amp, make_plot=make_plot)
            print("%i" % full_well)
        except RuntimeError:
            full_well = detResp.full_well(amp, frac_offset=0.05,
                                          make_plot=make_plot)
            print("%i (5%% dev)" % full_well)
