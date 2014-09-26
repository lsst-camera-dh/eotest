import os
import sys
import glob
from collections import OrderedDict
import numpy as np
import pyfits
import pylab_plotter as plot
from EOTestResults import EOTestResults
from fe55_gain_fitter import fe55_gain_fitter
from DetectorResponse import DetectorResponse
import lsst.eotest.image_utils as imutils

class EOTestPlots(object):
    plotter = plot
    def __init__(self, sensor_id, rootdir='.', output_dir='.', ps=False,
                 interactive=False):
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
        results_file = self._fullpath('%s_eotest_results.fits' % sensor_id)
        if not os.path.exists(results_file):
            raise RuntimeError("EOTestPlots: %s not found" % results_file)
        self.results = EOTestResults(results_file)
        self._qe_data = None
    @property
    def qe_data(self):
        if self._qe_data is None:
            self._qe_data = pyfits.open(self._fullpath('%s_QE.fits' 
                                                       % self.sensor_id))
        return self._qe_data
    def _save_fig(self, outfile_root):
        plot.pylab.savefig(self._outputpath('%s.png' % outfile_root))
        if self.ps:
            plot.pylab.savefig(self._outputpath('%s.eps' % outfile_root))
    def _fullpath(self, basename):
        return os.path.join(self.rootdir, basename)
    def _outputpath(self, basename):
        return os.path.join(self.output_dir, basename)
    def fe55_dists(self, chiprob_min=0.1):
        fe55_file = glob.glob(self._fullpath('%s_psf_results*.fits' 
                                                  % self.sensor_id))[0]
        fe55_catalog = pyfits.open(fe55_file)
        for amp in imutils.allAmps:
            print "Amp", amp
            chiprob = fe55_catalog[amp].data.field('CHIPROB')
            index = np.where(chiprob > chiprob_min)
            dn = fe55_catalog[amp].data.field('DN')[index]
            gain, peak, sigma = fe55_gain_fitter(dn, make_plot=True,
                                                 title='Amp %i' % amp,
                                                 interactive=self.interactive)
            self._save_fig('Fe55_dist_%s_amp%02i' % (self.sensor_id, amp))
    def ptcs(self, xrange=(0.1, 1e4), yrange=(0.1, 1e4)):
        ptc = pyfits.open(self._fullpath('%s_ptc.fits' % self.sensor_id))
        for amp in imutils.allAmps:
            print "Amp", amp
            mean = ptc[1].data.field('AMP%02i_MEAN' % amp)
            var = ptc[1].data.field('AMP%02i_VAR' % amp)
            win = plot.xyplot(mean, var, xname='mean (ADU)',
                              yname='variance (ADU^2)',
                              xrange=xrange, yrange=yrange,
                              xlog=1, ylog=1)
            xx = np.logspace(np.log10(xrange[0]), np.log10(xrange[1]), 20)
            plot.curve(xx, xx/self.results['GAIN'][amp-1], oplot=1, color='r')
            win.set_title('Amp %i' % amp)
            self._save_fig('PTC_%s_amp%02i' % (self.sensor_id, amp))
    def gains(self):
        results = self.results
        plot.xyplot(results['AMP'], results['GAIN'], xname='Amp',
                    yname='gain (e-/ADU)', yrange=(0, max(results['GAIN'])*1.2))
        self._save_fig('Gain_vs_amp_%s' % self.sensor_id)
    def linearity(self, gain_range=(1, 6), max_dev=0.02):
        ptc = pyfits.open(self._fullpath('%s_ptc.fits' % self.sensor_id))
        detresp = DetectorResponse(self._fullpath('%s_det_response.fits' 
                                                  % self.sensor_id),
                                   ptc=ptc, gain_range=gain_range)
        for amp in imutils.allAmps:
            print "Amp", amp
            try:
                maxdev, fit_pars = detresp.linearity(amp, make_plot=True,
                                                     title='Amp %i' % amp,
                                                     max_dev=max_dev,
                                                     interactive=self.interactive)
                self._save_fig('Linearity_%s_amp%02i' % (self.sensor_id, amp))
            except Exception, e:
                print e
                continue
    def qe_ratio(self, ref, amp=None):
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
            self._save_fig('QE_ratio_%s_amp%02i' % (self.sensor_id, amp))
    def qe(self):
        qe_data = self.qe_data
        band_pass = OrderedDict()
        band_pass['u'] = (321, 391)
        band_pass['g'] = (402, 552)
        band_pass['r'] = (552, 691)
        band_pass['i'] = (691, 818)
        band_pass['z'] = (818, 922)
        band_pass['y'] = (930, 1070)
        band_wls = np.array([sum(band_pass[b])/2. for b in band_pass.keys()])

        wl = qe_data[1].data.field('WAVELENGTH')
        qe = {}
        for amp in imutils.allAmps:
            qe[amp] = qe_data[1].data.field('AMP%02i' % amp)
            win = plot.curve(wl, qe[amp], xname='wavelength (nm)',
                             yname='QE (e-/photon)', oplot=amp-1)
            qe_band = qe_data[2].data.field('AMP%02i' % amp)
            plot.xyplot(band_wls, qe_band, oplot=1)
        self._save_fig('QE_%s' % self.sensor_id)
    def confluence_table(self, outfile=False):
        if outfile:
            output = open(self._outputpath('%s_results.txt' % self.sensor_id), 'w')
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

if __name__ == '__main__':
    plots = EOTestPlots('114-03')
#    plots.fe55_dists(chiprob_min=0.1)
#    plots.ptcs()
#    plots.gains()
#    plots.linearity()
    plots.confluence_table()
