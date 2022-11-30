"""
Find the footprints of Fe55 clusters and do statistics on the pixel
values.
"""
from __future__ import absolute_import, print_function, division
import pickle
import numpy as np
import numpy.lib.recfunctions as nlr
import matplotlib.pyplot as plt
import scipy.stats
import scipy.optimize
import lsst.afw.detection as afwDetect
import lsst.afw.math as afwMath
from .MaskedCCD import MaskedCCD
from .fe55_gain_fitter import fe55_lines

__all__ = ['Fe55PixelStats']


class RecArray(object):
    """
    Wrapper class to give numpy.recarray custom formatted __repr__ output.
    """
    _formats = dict((('<i4', '%12i'), ('<i8', '%12i'), ('<f8', '%12.3e')))

    def __init__(self, recarray):
        "Constructor"
        self.recarray = recarray

    def __repr__(self):
        """
        Formatting with column name.
        """
        lines = [''.join(['%12s' % name for name in self.dtype.names])]
        columns = tuple([self[colname] for colname in self.dtype.names])
        dtypes = tuple(self._formats[x[1]] for x in self.dtype.descr)
        for row in zip(*columns):
            lines.append(''.join([fmt % item for fmt, item
                                  in zip(dtypes, row)]))
        return '\n'.join(lines)

    def __getitem__(self, key):
        "Access column by key."
        return self.recarray[key]

    def __getattr__(self, attr):
        "Delegate everything else to self.recarray"
        return getattr(self.recarray, attr)

    @staticmethod
    def create(data, names=None):
        "Factory method wrapping numpy.rec.array"
        return RecArray(np.rec.array(data, names=names))


def profile_plot(ax, xarg, yarg, bins=20, xbounds=None, color='black',
                 plot_points=True):
    """
    Create a profile plot and return a numpy recarray of the plotted
    profile.
    """
    if xbounds is None:
        xmin, xmax = min(xarg), max(xarg)
    else:
        xmin, xmax = xbounds
    bin_edges = np.linspace(xmin, xmax, bins+1)
    xx = np.array(xarg)
    yy = np.array(yarg)
    bin_indices = np.digitize(xx, bin_edges)
    data = []
    for bin_ in range(bins):
        index = np.where(bin_indices == bin_)
        if len(index[0]) == 0:
            continue
        data.append(((bin_edges[bin_] + bin_edges[bin_+1])/2.,
                     np.median(yy[index]), np.std(yy[index])))
    my_recarr = np.rec.array(data, names='bin_centers ymedian yerr'.split())
    if plot_points:
        ax.errorbar(my_recarr['bin_centers'], my_recarr['ymedian'],
                    yerr=my_recarr['yerr'], xerr=(xmax-xmin)/float(bins)/2.,
                    fmt="none", ecolor=color)
    return my_recarr


def get_fp_pixels(ccd, amp, nsig=4, bg_reg=(10, 10), npix_range=(5, 20)):
    """
    Return a numpy record array with rows of the 9 pixel values around
    the peaks of Fe55 clusters.  The footprint threshold is
    nsig*stdevclip + median, bg_reg is the NxM local background
    region, and the number of pixels per footprint is restricted to be
    within npix_range.
    """
    # Background subtraction using a background image based on the
    # local background levels.
    bg_ctrl = afwMath.BackgroundControl(bg_reg[0], bg_reg[1], ccd.stat_ctrl)
    image = ccd[amp]
    image -= afwMath.makeBackground(ccd[amp], bg_ctrl).getImageF()

    # Set the footprint threshold.
    stats = afwMath.makeStatistics(image, afwMath.MEDIAN | afwMath.STDEVCLIP,
                                   ccd.stat_ctrl)
    threshold = afwDetect.Threshold(nsig*stats.getValue(afwMath.STDEVCLIP)
                                    + stats.getValue(afwMath.MEDIAN))

    # Gather the data for the record array object.  This includes the
    # amplifier number, the coordinates of the peak pixel (x0, y0),
    # the 3x3 pixels centered on (x0, y0), and the sum of ADU values
    # over the footprint.
    imarr = image.getImage().getArray()
    data = []
    for fp in afwDetect.FootprintSet(image, threshold).getFootprints():
        if fp.getArea() < npix_range[0] or npix_range[1] < fp.getArea():
            continue
        # Get peak coordinates and p0-p8 values
        peaks = [pk for pk in fp.getPeaks()]
        if len(peaks) > 1:
            continue
        x0, y0 = peaks[0].getIx(), peaks[0].getIy()
        row = [amp, x0, y0] + list(imarr[y0-1:y0+2, x0-1:x0+2].flatten())
        if len(row) != 12:
            # Skip if there are too few p0-p8 values (e.g., cluster is
            # on the imaging section edge).
            continue
        # Sum over pixels in the footprint.
        dn_sum = 0
        spans = fp.getSpans()
        for span in spans:
            y = span.getY()
            for x in range(span.getX0(), span.getX1() + 1):
                dn_sum += imarr[y][x]
        row.append(dn_sum)
        data.append(row)
    names = 'amp x y'.split() + ['p%i' % i for i in range(9)] + ['DN_sum']

    return np.rec.array(data, names=names)


class Fe55PixelStats(object):
    "Statistics of pixel values around an Fe55 cluster peak."
    _selections = dict(amp='self._amp_selection',
                       kalpha='self._kalpha_selection')

    def __init__(self, input_files, mask_files=(), sensor_id=None,
                 logger=None, selection='amp', linearity_correction=None):
        """
        Extract record array from the record arrays for all of the
        input_files.
        """
        self.sensor_id = sensor_id
        self.set_selection_function(selection)
        rec_arrays = []
        for infile in input_files:
            if logger is not None:
                logger.info("Processing %s" % infile)
            ccd = MaskedCCD(infile, mask_files=mask_files,
                            linearity_correction=linearity_correction)
            for amp in ccd:
                rec_arrays.append(get_fp_pixels(ccd, amp))
        self.rec_array = nlr.stack_arrays(rec_arrays, usemask=False,
                                          autoconvert=True, asrecarray=True)
        self.amps = sorted(ccd.keys())

    def set_selection_function(self, selection):
        if selection not in self._selections:
            raise RuntimeError("Unrecognized data selection: " + selection)
        self._selection_key = selection

    @property
    def _selection(self):
        return eval(self._selections[self._selection_key])

    def to_pickle(self, outfile):
        pickle.dump(self, open(outfile, 'wb'))

    @staticmethod
    def read_pickle(infile):
        return pickle.load(open(infile, 'rb'))

    def _multi_panel_figure(self, figsize):
        """
        Boilerplate code for setting up the figure for showing plots
        for all of the amplifiers.
        """
        plt.rcParams['figure.figsize'] = figsize
        fig = plt.figure()
        if self.sensor_id is not None:
            frame_axes = fig.add_subplot(111, frameon=False)
            frame_axes.set_title(self.sensor_id)
        return fig, frame_axes

    def _amp_selection(self, amp, **kwds):
        """
        Select rows based only on amplifier.
        """
        return np.where(self.rec_array.amp == amp)

    def _kalpha_selection(self, amp, bins=50, nsig=2):
        """
        For the specified amplifier, fit the DN distribution to
        idenitfy the clusters in the K-alpha peak.
        """
        my_recarr = self.rec_array[np.where(self.rec_array.amp == amp)]
        dn = my_recarr['DN_sum']
        median = np.median(dn)
        #stdev = np.std(dn)
        stdev = afwMath.makeStatistics(np.array(dn, dtype=float),
                                       afwMath.STDEVCLIP).getValue()
        dn_range = (median - nsig*stdev, median + nsig*stdev)
        results = np.histogram(dn, bins=bins, range=dn_range)
        x = (results[1][1:] + results[1][:-1])/2.
        y = results[0]
        ntot = sum(y)
        p0 = (0.88*ntot, median, stdev/2., 0.12*ntot)
        pars, _ = scipy.optimize.curve_fit(fe55_lines, x, y, p0=p0)
        dn_min, dn_max = pars[1] - nsig*pars[2], pars[1] + nsig*pars[2]
        index = np.where((self.rec_array.amp == amp) &
                         (self.rec_array.DN_sum > dn_min) &
                         (self.rec_array.DN_sum < dn_max))
        return index

    @staticmethod
    def apply_offsets(axes, xoffset=0.025, yoffset=0.025):
        bbox = axes.get_position()
        points = bbox.get_points()
        points[0] += xoffset
        points[1] += yoffset
        bbox.set_points(points)
        axes.set_position(bbox)

    def pixel_hists(self, pix0='p3', pix1='p5', figsize=(10, 10), bins=50,
                    dn_range=(-10, 30)):
        """
        Plot histograms of pix0 and pix1 values.
        """
        fig, frame_axes = self._multi_panel_figure(figsize)
        frame_axes.set_xlabel('%(pix0)s (blue), %(pix1)s (red) (ADU)'
                              % locals())
        frame_axes.get_xaxis().set_ticks([])
        frame_axes.get_yaxis().set_ticks([])
        for amp in self.amps:
            subplot = (4, 4, amp)
            ax = fig.add_subplot(*subplot)
            self.apply_offsets(ax)
            selection = self._selection(amp, bins=bins)
            my_recarr = self.rec_array[selection]
            plt.hist(my_recarr[pix0], color='blue', histtype='step',
                     range=dn_range, bins=20)
            plt.hist(my_recarr[pix1], color='red', histtype='step',
                     range=dn_range, bins=20)
            if amp in (1, 5, 9, 13):
                ax.set_ylabel('entries / bin')
            plt.annotate('Amp %i' % amp, (0.5, 0.9),
                         xycoords='axes fraction', size='x-small')
        return fig

    def pixel_diff_profile(self, pixel_coord='x', pix0='p3', pix1='p5',
                           bins=50, figsize=(10, 10)):
        """
        Fit the difference of pix1 and pix0 profiles as a function of
        cluster peak pixel coordinate.
        """
        fig, frame_axes = self._multi_panel_figure(figsize)
        frame_axes.set_xlabel('%s pixel index (%s, red; %s, blue)'
                              % (pixel_coord, pix1, pix0))
        frame_axes.get_xaxis().set_ticks([])
        frame_axes.get_yaxis().set_ticks([])
        data = []
        for amp in self.amps:
            subplot = (4, 4, amp)
            axes = fig.add_subplot(*subplot)
            self.apply_offsets(axes)
            selection = self._selection(amp, bins=bins)
            my_recarr = self.rec_array[selection]
            p1_prof = profile_plot(axes, my_recarr[pixel_coord],
                                   my_recarr[pix1], color='red',
                                   plot_points=True)
            p0_prof = profile_plot(axes, my_recarr[pixel_coord],
                                   my_recarr[pix0], color='blue',
                                   plot_points=True)
            plt.annotate('Amp %i' % amp, (0.5, 0.9),
                         xycoords='axes fraction', size='x-small')
            plt.xlim(min(self.rec_array[pixel_coord]),
                     max(self.rec_array[pixel_coord]))
            # Fit a line to the pix1 + pix0 profile.
            x = p1_prof['bin_centers']
            y = (p1_prof['ymedian'] + p0_prof['ymedian'])/2.
            pars, cov = np.polyfit(x, y, 1, cov=True)
            error = np.sqrt(cov[0][0])
            ymodel = np.poly1d(pars)(x)
            plt.plot(x, ymodel, color='green')
            if amp in (1, 5, 9, 13):
                axes.set_ylabel('DN (bg-subtracted)')
            data.append([amp, pars[0], error, pars[0]/error])
        return fig, RecArray.create(data, names='amp slope error nsig'.split())

    def apflux_profile(self, pixel_coord='x', bins=50, figsize=(10, 10)):
        """
        Plot the flux summed over the 9 pixels as a function of pixel
        coordinate.
        """
        fig, frame_axes = self._multi_panel_figure(figsize)
        frame_axes.set_xlabel('%s pixel index' % pixel_coord)
        frame_axes.get_xaxis().set_ticks([])
        frame_axes.get_yaxis().set_ticks([])
        data = []
        for amp in self.amps:
            subplot = (4, 4, amp)
            axes = fig.add_subplot(*subplot)
            self.apply_offsets(axes)
            selection = self._selection(amp, bins=bins)
            my_recarr = self.rec_array[selection]
            apflux = sum(my_recarr['p%i' % i] for i in range(9))
            apflux_prof = profile_plot(axes, my_recarr[pixel_coord],
                                       apflux, plot_points=True)
            plt.annotate('Amp %i' % amp, (0.5, 0.9),
                         xycoords='axes fraction', size='x-small')
            plt.xlim(min(self.rec_array[pixel_coord]),
                     max(self.rec_array[pixel_coord]))
            # Fit a line to the apflux profile.
            x = apflux_prof['bin_centers']
            y = apflux_prof['ymedian']
            pars, cov = np.polyfit(x, y, 1, cov=True)
            error = np.sqrt(cov[0][0])
            func = np.poly1d(pars)
            plt.plot(x, func(x), color='green')
            if amp in (1, 5, 9, 13):
                axes.set_ylabel('Aperture flux (ADU)')
            data.append([amp, pars[0], error, pars[0]/error])
        return fig, RecArray.create(data, names='amp slope error nsig'.split())

    def dn_hists(self, figsize=(10, 10), nsig=2, bins=30):
        """
        Plot histograms of the pixel DNs summed over each cluster
        footprint.
        """
        fig, frame_axes = self._multi_panel_figure(figsize)
        frame_axes.set_xlabel('Pixel DNs summed over cluster footprint')
        frame_axes.get_xaxis().set_ticks([])
        frame_axes.get_yaxis().set_ticks([])
        for amp in self.amps:
            subplot = (4, 4, amp)
            ax = fig.add_subplot(*subplot)
            self.apply_offsets(ax)
            my_recarr = self.rec_array[np.where(self.rec_array.amp == amp)]
            dn = my_recarr['DN_sum']
            median = np.median(dn)
            stdev = afwMath.makeStatistics(np.array(dn, dtype=float),
                                           afwMath.STDEVCLIP).getValue()
            dn_range = (median - nsig*stdev, median + nsig*stdev)
            plt.hist(dn, histtype='step', bins=bins, range=dn_range)
            selection = self._selection(amp, bins=bins)
            plt.hist(self.rec_array[selection]['DN_sum'], histtype='step',
                     bins=bins, range=dn_range, color='red')
            if amp in (1, 5, 9, 13):
                ax.set_ylabel('entries / bin')
            plt.annotate('Amp %i' % amp, (0.5, 0.9),
                         xycoords='axes fraction', size='x-small')
        return fig
