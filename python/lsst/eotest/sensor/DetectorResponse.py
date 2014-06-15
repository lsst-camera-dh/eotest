"""
@brief Class to compute detector response characteristics (e.g., full well,
linearity) from flat pairs data.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import sys
import numpy as np
import pyfits
import scipy.optimize
import lsst.eotest.image_utils as imutils
import pylab_plotter as plot

class DetectorResponse(object):
    def __init__(self, infile, amps=imutils.allAmps, ptc=None,
                 gain_range=None):
        if infile[-5:] == '.fits':
            self._read_from_fits(infile)
        else:
            self._read_from_text(infile)
        self._sort_by_flux()
        self._compute_gain_selection(ptc, gain_range)
    def _compute_gain_selection(self, ptc, gain_range):
        self._index = {}
        if ptc is None or gain_range is None:
            return
        if len(ptc[1].data.field('AMP01_MEAN')) != len(self.flux):
            raise RuntimeError('Number of measurements in PTC file ' +
                               'differs from detector response file.')
        for amp in imutils.allAmps:
            mean = ptc[1].data.field('AMP%02i_MEAN' % amp)
            var = ptc[1].data.field('AMP%02i_VAR' % amp)
            gain = mean/var
            self._index[amp] = np.where((gain >= gain_range[0]) & 
                                        (gain <= gain_range[1]))
    def _sort_by_flux(self):
        index = np.argsort(self.flux)
        self.flux = self.flux[index]
        for amp in imutils.allAmps:
            self.Ne[amp] = self.Ne[amp][index]
    def _read_from_fits(self, infile):
        foo = pyfits.open(infile)
        hdu = foo['DETECTOR_RESPONSE']
        self.flux = np.array(hdu.data.field('FLUX'), dtype=np.float)
        self.Ne = dict([(amp, np.array(hdu.data.field('AMP%02i_SIGNAL' % amp),
                                       dtype=np.float))
                        for amp in imutils.allAmps])
    def _read_from_text(self, infile):
        data = np.recfromtxt(infile)
        data = data.transpose()
        self.flux = data[0]
        self.Ne = dict([(amp, ne) for amp, ne in
                        zip(imutils.allAmps, data[1:])])
    def full_well(self, amp, order=15, fit_range=(1e2, 5e4),
                  poly_fit_min=1e3, frac_offset=0.1, dfrac=None,
                  make_plot=False):
        flux = self.flux
        Ne = self.Ne[amp]
        if self._index:
            flux = flux[self._index[amp]]
            Ne = Ne[self._index[amp]]
        #
        # Fit linear part of response curve.
        #
        indx1 = np.where((Ne > fit_range[0]) & (Ne < fit_range[1]))
        f1 = np.poly1d(np.polyfit(flux[indx1], Ne[indx1], 1))
        #
        # Fit polynomial of specified order.
        #
        Ne_vals = f1(flux)
        dNe_frac = (Ne - Ne_vals)/Ne_vals
        if dNe_frac[-1] > -frac_offset:
            message = """Detector response data do not extend past the
              specified fractional offset of %i%% for the full
              well definition.""" % (int(frac_offset*100),)
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
        df = lambda xx : 1 - fp(xx)/f1(xx) - frac_offset
        x = flux[indxp]
        flux0 = scipy.optimize.brentq(df, x[len(x)/2], x[-1])
        full_well = int(fp(flux0))
        
        # Save pylab interactive state.
        pylab_interactive_state = plot.pylab.isinteractive()
        if make_plot:
            plot.pylab.ion()
            yrange = (0, max(max(Ne), full_well)*1.1)
            plot.xyplot(flux, Ne,
                        xname='Illumination (flux*exptime, arb. units)',
                        yname='e- per pixel',
                        yrange=yrange)
            plot.curve(flux, f1(flux), oplot=1, color='g')
            plot.xyplot(flux[indxp], Ne[indxp], oplot=1, color='r')
            xx = np.linspace(0, max(flux), 100)
            plot.curve(xx, fp(xx), oplot=1, color='b')
            plot.vline(flux0)
            plot.hline(full_well)
            plot.pylab.annotate('Segment %s\nfull well = %i'
                                % (imutils.channelIds[amp], int(full_well)),
                                (0.1, 0.8), xycoords='axes fraction')
        # Restore pylab interactive state.
        plot.pylab.interactive(pylab_interactive_state)
        return full_well, fp
    def linearity(self, amp, fit_range=(1e2, 9e4), max_dev=0.02,
                  make_plot=False, title=None, interactive=True):
        flux, Ne = self.flux, self.Ne[amp]
        if self._index:
            flux = flux[self._index[amp]]
            Ne = Ne[self._index[amp]]
        indx = np.where((Ne > fit_range[0]) & (Ne < fit_range[1]))
        f1_pars = np.polyfit(flux[indx], Ne[indx], 1)
        f1 = np.poly1d(f1_pars)
        dNfrac = 1 - Ne/f1(flux)
        if make_plot:
            if interactive:
                plot.pylab.ion()
            else:
                plot.pylab.ioff()
            win = plot.xyplot(Ne, dNfrac, xname='e-/pixel', xlog=1,
                              yname='frac. deviation of e-/pixel from linear fit',
                              xrange=(1e2, 2e5),
                              yrange=(-1.5*max_dev, 1.5*max_dev))
            plot.curve(Ne, dNfrac, oplot=1, color='r')
            plot.pylab.annotate("slope = %.2e\n\nintercept = %.2e" 
                                % tuple(f1_pars),
                                (0.1, 0.8), xycoords='axes fraction')
            if title is not None:
                win.set_title(title)
            plot.hline(0, lineStyle='--')
            plot.hline(max_dev)
            plot.hline(-max_dev)
            plot.vline(Ne[indx][0])
            plot.vline(Ne[indx][-1])
        return max(abs(dNfrac[indx])), f1_pars

if __name__ == '__main__':
    infile = '000-00_det_response.txt'
    detResp = DetectorResponse(infile)
    make_plot = False
    print "Segment   max. dev.   full well"
    for amp in imutils.allAmps:
        max_dev, fit_pars = detResp.linearity(amp, make_plot=make_plot)
        sys.stdout.write("%s         %.3f       "
                         % (imutils.channelIds[amp], max_dev))
        try:
            full_well =  detResp.full_well(amp, make_plot=make_plot)
            print "%i" % full_well
        except RuntimeError:
            full_well = detResp.full_well(amp, frac_offset=0.05,
                                          make_plot=make_plot)
            print "%i (5%% dev)" % full_well
