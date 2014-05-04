import os
import sys
import numpy as np
import pyfits
import pylab_plotter as plot
import lsst.eotest.sensor as sensorTest
import lsst.eotest.image_utils as imutils

plot.pylab.ion()

class EOTestPlots(object):
    def __init__(self, sensor_id):
        self.sensor_id = sensor_id
        self.results = sensorTest.EOTestResults('%s_eotest_results.fits' 
                                                % sensor_id)
    def fe55_dists(self, chiprob_min=0.1):
        fe55_catalog = pyfits.open('%s_psf_results.fits' % self.sensor_id)
        for amp in imutils.allAmps:
            chiprob = fe55_catalog[amp].data.field('CHIPROB')
            index = np.where(chiprob > chiprob_min)
            dn = fe55_catalog[amp].data.field('DN')[index]
            gain, peak, sig = sensorTest.fe55_gain_fitter(dn, make_plot=True,
                                                          title='Amp %i' % amp)
        plot.pylab.savefig('Fe55_dist_%s_amp%02i.png' % (self.sensor_id, amp))
    def ptcs(self, xrange=(0.1, 1e4), yrange=(0.1, 1e4)):
        ptc = pyfits.open('%s_ptc.fits' % self.sensor_id)
        for amp in imutils.allAmps:
            mean = ptc[1].data.field('AMP%02i_MEAN' % amp)
            var = ptc[1].data.field('AMP%02i_VAR' % amp)
            win = plot.xyplot(mean, var, xname='mean (ADU)',
                              yname='variance (ADU^2)',
                              xrange=xrange, yrange=yrange,
                              xlog=1, ylog=1)
        xx = np.logspace(np.log10(xrange[0]), np.log10(xrange[1]), 20)
        plot.curve(xx, xx/self.results['GAIN'][amp-1], oplot=1, color='r')
        win.set_title('Amp %i' % amp)
        plot.pylab.savefig('PTC_%s_amp%02i.png' % (self.sensor_id, amp))
    def gains(self):
        plot.xyplot(results['AMP'], results['GAIN'], xname='Amp',
                    yname='gain (e-/ADU)', yrange=(0, max(results['GAIN'])*1.2))
        plot.pylab.savefig('Gain_vs_amp_%s.png' % self.sensor_id)
    def linearity(self, gain_range=(1, 6), max_dev=0.02):
        ptc = pyfits.open('%s_ptc.fits' % self.sensor_id)
        detresp = sensorTest.DetectorResponse('%s_det_response.fits' 
                                              % self.sensor_id,
                                              ptc=ptc, gain_range=gain_range)
        for amp in imutils.allAmps:
            try:
                maxdev, fit_pars = detresp.linearity(amp, make_plot=True,
                                                     title='Amp %i' % amp,
                                                     max_dev=max_dev)
                plot.pylab.savefig('Linearity_%s_amp%02i.png' 
                                   % (self.sensor_id, amp))
            except Exception, e:
                print e
                continue
    def confluence_table(self, outfile=False):
        if outfile:
            output = open('%s_results.txt' % self.sensor_id, 'w')
        else:
            output = sys.stdout
        for name in self.results.colnames:
            output.write('|| %s' % name)
        output.write('||\n')
        for i, amp in enumerate(results['AMP']):
            output.write('| %i | %.2f | %.2f | %i | %.1e | %.1e | %.1e | %i | %i |\n'
                         % tuple([results[x][i] for x in results.colnames]))
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
