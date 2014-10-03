import os
import sys
import glob
import copy
from collections import OrderedDict
import numpy as np
import pyfits
import pylab
import pylab_plotter as plot
from EOTestResults import EOTestResults
from Fe55GainFitter import Fe55GainFitter
from DetectorResponse import DetectorResponse
from crosstalk import CrosstalkMatrix
import lsst.eotest.image_utils as imutils

class EOTestPlots(object):
    plotter = plot
    def __init__(self, sensor_id, rootdir='.', output_dir='.', ps=False,
                 interactive=False, results_file=None):
        self.sensor_id = sensor_id
        self.rootdir = rootdir
        self.output_dir = output_dir
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        self.ps = ps
        self.interactive = interactive
        if interactive:
            plot.pylab.ion()
        else:
            plot.pylab.ioff()
        if results_file is None:
            results_file = self._fullpath('%s_eotest_results.fits' % sensor_id)
        if not os.path.exists(results_file):
            raise RuntimeError("EOTestPlots: %s not found" % results_file)
        self.results = EOTestResults(results_file)
        self._qe_data = None
        self._qe_file = self._fullpath('%s_QE.fits' % self.sensor_id)
    @property
    def qe_data(self):
        if self._qe_data is None:
            self._qe_data = pyfits.open(self._qe_file)
        return self._qe_data
    def _save_fig(self, outfile_root):
        plot.pylab.savefig(self._outputpath('%s.png' % outfile_root))
        if self.ps:
            plot.pylab.savefig(self._outputpath('%s.eps' % outfile_root))
    def _fullpath(self, basename):
        return os.path.join(self.rootdir, basename)
    def _outputpath(self, basename):
        return os.path.join(self.output_dir, basename)
    def crosstalk_matrix(self, cmap=pylab.cm.hot, xtalk_file=None):
        if xtalk_file is None:
            xtalk_file = os.path.join(self.rootdir, 
                                      '%s_xtalk_matrix.fits' % self.sensor_id)
        foo = CrosstalkMatrix(xtalk_file)
#        foo.plot_matrix(cmap=cmap)
        win = foo.plot(title="Crosstalk, %s" % self.sensor_id)
        return foo
    def fe55_dists(self, chiprob_min=0.1, fe55_file=None):
        if fe55_file is None:
            fe55_file = glob.glob(self._fullpath('%s_psf_results*.fits' 
                                                 % self.sensor_id))[0]
        fe55_catalog = pyfits.open(fe55_file)
        win = None
        figsize = (11, 8.5)
        for amp in imutils.allAmps:
            #print "Amp", amp
            chiprob = fe55_catalog[amp].data.field('CHIPROB')
            index = np.where(chiprob > chiprob_min)
            dn = fe55_catalog[amp].data.field('DN')[index]
            foo = Fe55GainFitter(dn)
            foo.fit()
            if amp == 1:
                win = foo.plot(interactive=True, subplot=(4, 4, amp),
                               figsize=figsize, frameLabels=True, amp=amp)
                win.frameAxes.text(0.5, 1.08, 'Fe55, %s' % self.sensor_id,
                                   horizontalalignment='center',
                                   verticalalignment='top',
                                   transform=win.frameAxes.transAxes,
                                   size='large')
            else:
                foo.plot(interactive=True, subplot=(4, 4, amp), win=win,
                         frameLabels=True, amp=amp)
            pylab.locator_params(axis='x', nbins=4, tight=True)
    def ptcs(self, xrange=(0.1, 1e4), yrange=(0.1, 1e4), figsize=(11, 8.5),
             ptc_file=None):
        if ptc_file is not None:
            ptc = pyfits.open(ptc_file)
        else:
            ptc = pyfits.open(self._fullpath('%s_ptc.fits' % self.sensor_id))
        for amp in imutils.allAmps:
            #print "Amp", amp
            subplot = (4, 4, amp)
            if amp == 1:
                win = plot.Window(subplot=subplot, figsize=figsize,
                                  xlabel=r'mean (ADU)',
                                  ylabel=r'variance (ADU$^2$)', size='large')
                win.frameAxes.text(0.5, 1.08,
                                   'Photon Transfer Curves, %s' \
                                       % self.sensor_id,
                                   horizontalalignment='center',
                                   verticalalignment='top',
                                   transform=win.frameAxes.transAxes,
                                   size='large')
            else:
                win.select_subplot(*subplot)
            self._offset_subplot(win)
            mean = ptc[1].data.field('AMP%02i_MEAN' % amp)
            var = ptc[1].data.field('AMP%02i_VAR' % amp)
            win = plot.xyplot(mean, var, xname='', yname='',
                              xrange=xrange, yrange=yrange,
                              xlog=1, ylog=1, new_win=False,)
            xx = np.logspace(np.log10(xrange[0]), np.log10(xrange[1]), 20)
            plot.curve(xx, xx/self.results['GAIN'][amp-1], oplot=1, color='r')
            pylab.annotate('Amp %i' % amp, (0.475, 0.8),
                           xycoords='axes fraction', size='x-small')
    def _offset_subplot(self, win, xoffset=0.025, yoffset=0.025):
        bbox = win.axes[-1].get_position()
        points = bbox.get_points()
        points[0] += xoffset
        points[1] += yoffset
        bbox.set_points(points)
        win.axes[-1].set_position(bbox)
    def gains(self, oplot=0, xoffset=0.25, width=0.5, color='b'):
        results = self.results
        win = plot.bar(results['AMP'] - xoffset, results['GAIN'],
                       xname='Amp', yname='gain (e-/DN)',
                       yrange=(0, max(results['GAIN']*1.2)),
                       xrange=(0, 17), color=color, width=width)
        win.set_title(self.sensor_id)
    def noise(self, oplot=0, xoffset=0.25, width=0.5, color='b'):
        results = self.results
        win = plot.bar(results['AMP'] - xoffset, results['READ_NOISE'],
                       xname='Amp', yname='read noise (rms e-)',
                       yrange=(0, max(results['GAIN']*1.2)),
                       xrange=(0, 17), color=color, width=width)
        win.set_title(self.sensor_id)
    def linearity(self, gain_range=(1, 6), max_dev=0.02, figsize=(11, 8.5),
                  ptc_file=None, detresp_file=None):
        if ptc_file is not None:
            ptc = pyfits.open(ptc_file)
        else:
            ptc = pyfits.open(self._fullpath('%s_ptc.fits' % self.sensor_id))
        if detresp_file is not None:
            detresp = DetectorResponse(detresp_file, ptc=ptc,
                                       gain_range=gain_range)
        else:
            detresp = DetectorResponse(self._fullpath('%s_det_response.fits' 
                                                      % self.sensor_id),
                                       ptc=ptc, gain_range=gain_range)
        for amp in imutils.allAmps:
            #print "Amp", amp
            maxdev, fit_pars, Ne, flux = detresp.linearity(amp)
            f1 = np.poly1d(fit_pars)
            dNfrac = 1 - Ne/f1(flux)

            subplot = (4, 4, amp)
            if amp == 1:
                win = plot.Window(subplot=subplot, figsize=figsize,
                                  xlabel='e-/pixel', ylabel='', size='large')
                win.frameAxes.text(0.5, 1.08, 'Linearity, %s' % self.sensor_id,
                                   horizontalalignment='center',
                                   verticalalignment='top',
                                   transform=win.frameAxes.transAxes,
                                   size='large')
            else:
                win.select_subplot(*subplot)
            self._offset_subplot(win)
            # Resize subplot for plotting flux vs exposure.
            bbox = win.axes[-1].get_position()
            top_pts = bbox.get_points()
            bot_pts = copy.deepcopy(top_pts)
            dx, dy = top_pts[1] - top_pts[0]
            top_pts[0][1] += dy/4
            bbox.set_points(top_pts)
            win.axes[-1].set_position(bbox)
            win.axes[-1].loglog(Ne, flux, 'ko')
            win.axes[-1].loglog(f1(flux), flux, 'r-')
            pylab.annotate('Amp %i' % amp, (0.2, 0.8),
                           xycoords='axes fraction', size='x-small')
            if amp in (1, 5, 9, 13):
                win.axes[-1].set_ylabel('flux')
            for label in win.axes[-1].get_xticklabels():
                label.set_visible(False)

            # Add fractional residuals sub-subplot.
            bot_rect = [bot_pts[0][0], bot_pts[0][1], dx, dy/4.]
            bot_ax = win.fig.add_axes(bot_rect, sharex=win.axes[-1])
            bot_ax.semilogx(Ne, dNfrac, 'ko')
            bot_ax.semilogx(Ne, np.zeros(len(Ne)), 'r-')
            pylab.locator_params(axis='y', nbins=5, tight=True)
            plot.setAxis(yrange=(-1.5*max_dev, 1.5*max_dev))
    def qe_ratio(self, ref, amp=None, qe_file=None):
        if qe_file is not None:
            self._qe_file = qe_file
        if amp is None:
            amps = imutils.allAmps
        else:
            amps = (amp,)
        for amp in amps:
            print "Amp", amp
            wls = []
            ref_wls = ref.qe_data[1].data.field('WAVELENGTH')
            fluxes, ref_fluxes = [], []
            column = 'AMP%02i' % amp
            for i, wl in enumerate(self.qe_data[1].data.field('WAVELENGTH')):
                if wl in ref_wls and wl not in wls:
                    wls.append(wl)
                    fluxes.append(self.qe_data[1].data.field(column)[i])
            for i, wl in enumerate(ref_wls):
                if wl in wls:
                    ref_fluxes.append(ref.qe_data[1].data.field(column)[i])
            fluxes = np.array(fluxes)
            ref_fluxes = np.array(ref_fluxes)
            if self.interactive:
                plot.pylab.ion()
            win = plot.xyplot(wls, fluxes/ref_fluxes, xname='wavelength (nm)',
                              yname='QE(%s) / QE(%s)' % (self.sensor_id,
                                                         ref.sensor_id))
            win.set_title('Amp %i' % amp)
            plot.hline(1)
    def qe(self, qe_file=None):
        if qe_file is not None:
            self._qe_file = qe_file
        qe_data = self.qe_data
        bands = qe_data[2].data.field('BAND')

        band_pass = OrderedDict()
        band_pass['u'] = (321, 391)
        band_pass['g'] = (402, 552)
        band_pass['r'] = (552, 691)
        band_pass['i'] = (691, 818)
        band_pass['z'] = (818, 922)
        band_pass['y'] = (930, 1070)
        band_wls = np.array([sum(band_pass[b])/2. for b in band_pass.keys()
                             if b in bands])
        band_wls_errs =  np.array([(band_pass[b][1] - band_pass[b][0])/2. 
                                   for b in band_pass.keys() if b in bands])
        wl = qe_data[1].data.field('WAVELENGTH')
        qe = {}
        for amp in imutils.allAmps:
            qe[amp] = qe_data[1].data.field('AMP%02i' % amp)
            win = plot.curve(wl, qe[amp], xname='wavelength (nm)',
                             yname='QE (e-/photon)', oplot=amp-1,
                             xrange=(300, 1100))
            if amp == 1:
                win.set_title(self.sensor_id)
            qe_band = qe_data[2].data.field('AMP%02i' % amp)
            plot.xyplot(band_wls, qe_band, xerr=band_wls_errs, 
                        oplot=1, color='g')
    def confluence_table(self, outfile=False):
        if outfile:
            output = open(self._outputpath('%s_results.txt' 
                                           % self.sensor_id), 'w')
        else:
            output = sys.stdout
        for name in self.results.colnames:
            output.write('|| %s' % name)
        output.write('||\n')
        format = '| %i | %.2f | %.2f | %i | %.1e | %.1e | %.1e | %i | %i | %.1e |\n'
        for i, amp in enumerate(self.results['AMP']):
            output.write(format % tuple([self.results[x][i] 
                                         for x in self.results.colnames]))
        output.write('\n')
        if outfile:
            output.close()
    def latex_table(self, outfile):
        pass

if __name__ == '__main__':
    plots = EOTestPlots('114-03')
    plots.fe55_dists()
    plots.ptcs()
    plots.linearity()
    plots.gains()
    plots.noise()
    plots.qe()
    plots.crosstalk_matrix()
    plots.confluence_table()
    plots.latex_table()
