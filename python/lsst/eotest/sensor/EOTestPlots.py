import os
import sys
import glob
import copy
import numpy as np
import pyfits
import pylab
import matplotlib as mpl
import pylab_plotter as plot
from MaskedCCD import MaskedCCD
from EOTestResults import EOTestResults
from Fe55GainFitter import Fe55GainFitter
from DetectorResponse import DetectorResponse
from crosstalk import CrosstalkMatrix
from QE import QE_Data
from AmplifierGeometry import parse_geom_kwd
import lsst.eotest.image_utils as imutils
import lsst.afw.math as afwMath

def plot_flat(infile, nsig=3, cmap=pylab.cm.hot, win=None, subplot=(1, 1, 1),
              figsize=None, wl=None, gains=None):
    ccd = MaskedCCD(infile)
    foo = pyfits.open(infile)
    detsize = parse_geom_kwd(foo[1].header['DETSIZE'])
    nx = detsize['xmax']
    ny = detsize['ymax']
    mosaic = np.zeros((ny, nx), dtype=np.float)
    for ypos in range(2):
        for xpos in range(8):
            amp = ypos*8 + xpos + 1
            #
            # Compute subarray boundaries for this segment in the mosaicked
            # image array.
            datasec = parse_geom_kwd(foo[amp].header['DATASEC'])
            dx = np.abs(datasec['xmax'] - datasec['xmin']) + 1
            dy = np.abs(datasec['ymax'] - datasec['ymin']) + 1
            if ypos == 0:
                xmin = nx - (xpos + 1)*dx
                ymin = dy
            else:
                xmin = xpos*dx
                ymin = 0
            xmax = xmin + dx
            ymax = ymin + dy
            #
            # Extract the bias-subtracted image for this segment
            segment_image = ccd.unbiased_and_trimmed_image(amp)
            subarr = segment_image.getImage().getArray()
            #
            # Determine flips in x- and y-directions in order to
            # get the (1, 1) pixel in the lower right corner.
            detsec = parse_geom_kwd(foo[amp].header['DETSEC'])
            if detsec['xmax'] > detsec['xmin']:  # flip in x-direction
                subarr = subarr[:,::-1]
            if detsec['ymax'] > detsec['ymin']:  # flip in y-direction
                subarr = subarr[::-1,:]
            #
            # Convert from ADU to e-
            if gains is not None:
                subarr *= gains[amp-1]
            #
            # Set the subarray in the mosaicked image.
            mosaic[ymin:ymax, xmin:xmax] = subarr
    #
    # Set the color map to extend over the range median +/- stdev(clipped)
    # of the pixel values.
    pixel_data = mosaic.flatten()
    stats = afwMath.makeStatistics(pixel_data,
                                   afwMath.STDEVCLIP | afwMath.MEDIAN)
    median = stats.getValue(afwMath.MEDIAN)
    stdev = stats.getValue(afwMath.STDEVCLIP)
    vmin = max(min(pixel_data), median - nsig*stdev)
    vmax = min(max(pixel_data), median + nsig*stdev)
    if win is None:
        win = plot.Window(subplot=subplot, figsize=figsize,
                          xlabel='', ylabel='')
    else:
        win.select_subplot(*subplot)
    image = win.axes[-1].imshow(mosaic, interpolation='nearest', cmap=cmap)
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    image.set_norm(norm)
    if wl is None:
        # Extract wavelength from file
        wl = foo[0].header['MONOWL']
    win.axes[-1].set_title('%i nm' % wl)
    win.fig.colorbar(image)
    # Turn off tick labels for x- and y-axes
    pylab.setp(win.axes[-1].get_xticklabels(), visible=False)
    pylab.setp(win.axes[-1].get_yticklabels(), visible=False)
    return win

class EOTestPlots(object):
    band_pass = QE_Data.band_pass
    def __init__(self, sensor_id, rootdir='.', output_dir='.', ps=False,
                 interactive=False, results_file=None):
        self.sensor_id = sensor_id
        self.rootdir = rootdir
        self.output_dir = output_dir
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        self.ps = ps
        self.interactive = interactive
        plot.pylab.interactive(interactive)
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
    def fe55_dists(self, chiprob_min=0.1, fe55_file=None, figsize = (11, 8.5)):
        if fe55_file is None:
            fe55_file = glob.glob(self._fullpath('%s_psf_results*.fits' 
                                                 % self.sensor_id))[0]
        fe55_catalog = pyfits.open(fe55_file)
        win = None
        for amp in imutils.allAmps:
            #print "Amp", amp
            chiprob = fe55_catalog[amp].data.field('CHIPROB')
            index = np.where(chiprob > chiprob_min)
            dn = fe55_catalog[amp].data.field('DN')[index]
            foo = Fe55GainFitter(dn)
            try:
                foo.fit()
            except:
                continue
            if amp == 1:
                win = foo.plot(interactive=self.interactive, 
                               subplot=(4, 4, amp),
                               figsize=figsize, frameLabels=True, amp=amp)
                win.frameAxes.text(0.5, 1.08, 'Fe55, %s' % self.sensor_id,
                                   horizontalalignment='center',
                                   verticalalignment='top',
                                   transform=win.frameAxes.transAxes,
                                   size='large')
            else:
                foo.plot(interactive=self.interactive, 
                         subplot=(4, 4, amp), win=win,
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
        win.set_title("System Gain, %s" % self.sensor_id)
    def noise(self, oplot=0, xoffset=0.25, width=0.5, color='b'):
        results = self.results
        ymax = max(1.2*max(results['READ_NOISE']), 10)
        win = plot.bar(results['AMP'] - xoffset, results['READ_NOISE'],
                       xname='Amp', yname='read noise (rms e-)',
                       xrange=(0, 17), color=color, width=width,
                       yrange=(0, ymax))
        plot.hline(8)
        win.set_title("Read Noise, %s" % self.sensor_id)
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
            try:
                maxdev, fit_pars, Ne, flux = detresp.linearity(amp)
            except Exception, eObj:
                print "EOTestPlots.linearity: amp %i" % amp
                print "  ", eObj
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
            try:
                win.axes[-1].loglog(Ne, flux, 'ko')
            except Exception, eObj:
                print "EOTestPlots.linearity: amp %i" % amp
                print "  ", eObj
            try:
                win.axes[-1].loglog(f1(flux), flux, 'r-')
            except Exception, eObj:
                print "EOTestPlots.linearity: amp %i" % amp
                print "  ", eObj
            sys.stdout.flush()
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
        band_wls = np.array([sum(self.band_pass[b])/2. for b in 
                             self.band_pass.keys() if b in bands])
        band_wls_errs = np.array([(self.band_pass[b][1]-self.band_pass[b][0])/2.
                                  for b in self.band_pass.keys() if b in bands])
        wl = qe_data[1].data.field('WAVELENGTH')
        qe = {}
        for amp in imutils.allAmps:
            qe[amp] = qe_data[1].data.field('AMP%02i' % amp)
            win = plot.curve(wl, qe[amp], xname='wavelength (nm)',
                             yname='QE (% e-/photon)', oplot=amp-1,
                             xrange=(300, 1100))
            if amp == 1:
                win.set_title('QE, %s' % self.sensor_id)
            qe_band = qe_data[2].data.field('AMP%02i' % amp)
            plot.xyplot(band_wls, qe_band, xerr=band_wls_errs, 
                        oplot=1, color='g')
        plot.hline(100)
    def flat_fields(self, lambda_dir, nsig=3, cmap=pylab.cm.hot,
                    figsize=(11, 8.5)):
        glob_string = os.path.join(lambda_dir, '*_lambda_*.fits')
        print glob_string
        flats = sorted(glob.glob(glob_string))
        flats = [x for x in flats if x.find('bias') == -1]
        wls = []
        for flat in flats:
            wls.append(int(float(os.path.basename(flat).split('_')[2])))
        wls = np.array(wls)
        #print wls
        for i, band in enumerate(self.band_pass):
            mid_wl = sum(self.band_pass[band])/2.
            dwl = np.abs(wls - mid_wl)
            # Find the observed wavelength closest to the center of
            # the bandpass.
            target = np.where(dwl == min(dwl))[0][0]
            subplot = (2, 3, i+1)
            if i == 0:
                win = plot_flat(flats[target], subplot=subplot, nsig=nsig,
                                cmap=cmap, wl=wls[target], figsize=figsize,
                                gains=self.results['GAIN'])
            else:
                plot_flat(flats[target], subplot=subplot, win=win,
                          nsig=nsig, cmap=cmap, wl=wls[target],
                          gains=self.results['GAIN'])
        win.frameAxes.set_title('Flat Fields, %s' % self.sensor_id)
        return win
    def confluence_tables(self, outfile=False, prnu_file=None):
        if outfile:
            output = open(self._outputpath('%s_results.txt' 
                                           % self.sensor_id), 'w')
        else:
            output = sys.stdout
        # Write the per amp results.
        for name in self.results.colnames:
            output.write('|| %s' % name)
        output.write('||\n')
        format = '| %i | %.2f | %.2f | %i | %.1e | %.1e | %.1e | %i | %i | %.1e |\n'
        for i, amp in enumerate(self.results['AMP']):
            output.write(format % tuple([self.results[x][i] 
                                         for x in self.results.colnames]))
        output.write('\n')
        # Write the CCD-wide results.
        # PRNU:
        if prnu_file is None:
            prnu_file = self.results.infile
        prnu_results = pyfits.open(prnu_file)['PRNU_RESULTS'].data
        output.write("|| wavelength || stdev of pixel values || mean || stdev/mean ||\n")
        for wl, stdev, mean in zip(prnu_results['WAVELENGTH'],
                                   prnu_results['STDEV'], prnu_results['MEDIAN']):
            if stdev > 0:
                output.write("| %i | %12.4e | %12.4e | %12.4e |\n" 
                             % (wl, stdev, mean, stdev/mean))
            else:
                output.write("| %i | ... | ... | ... |\n" % wl)
        if outfile:
            output.close()
    def latex_tables(self, outfile):
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
