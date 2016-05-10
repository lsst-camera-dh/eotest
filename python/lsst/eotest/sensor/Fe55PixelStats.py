"""
Find the footprints of Fe55 clusters and do statistics on the pixel
values.
"""
from __future__ import absolute_import, print_function, division
import glob
import numpy as np
import numpy.lib.recfunctions as nlr
import matplotlib.pyplot as plt
import lsst.afw.detection as afwDetect
import lsst.afw.math as afwMath
from .MaskedCCD import MaskedCCD

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

def profile_plot(ax, xarg, yarg, bins=20, xbounds=None, color='black'):
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
                     np.mean(yy[index]), np.std(yy[index])))
    my_recarr = np.rec.array(data, names='bin_centers ymean yerr'.split())
    ax.errorbar(my_recarr['bin_centers'], my_recarr['ymean'],
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

    # Gather the data for the record array object.
    imarr = image.getImage().getArray()
    data = []
    for fp in afwDetect.FootprintSet(image, threshold).getFootprints():
        if fp.getNpix() < npix_range[0] or npix_range[1] < fp.getNpix():
            continue
        peaks = [pk for pk in fp.getPeaks()]
        if len(peaks) > 1:
            continue
        x0, y0 = peaks[0].getIx(), peaks[0].getIy()
        row = [amp, x0, y0] + list(imarr[y0-1:y0+2, x0-1:x0+2].flatten())
        if len(row) != 12:
            continue
        data.append(row)
    names = 'amp x y'.split() + ['p%i' % i for i in range(9)]
    return np.rec.array(data, names=names)

class Fe55PixelStats(object):
    "Statistics of pixel values around an Fe55 cluster peak."
    def __init__(self, input_files, mask_files=(), logger=None):
        """
        Extract record array from the record arrays for all of the
        input_files.
        """
        rec_arrays = []
        for infile in input_files:
            if logger is not None:
                logger.info("Processing %s" % infile)
            ccd = MaskedCCD(infile, mask_files=mask_files)
            for amp in ccd:
                rec_arrays.append(get_fp_pixels(ccd, amp))
        self.rec_array = nlr.stack_arrays(rec_arrays, usemask=False,
                                          autoconvert=True, asrecarray=True)
        self.amps = sorted(ccd.keys())

    def pixel_hists(self, pix0='p3', pix1='p5', figsize=(10, 10)):
        """
        Plot histograms of pix0 and pix1 values.
        """
        plt.rcParams['figure.figsize'] = figsize
        fig = plt.figure()
        for amp in self.amps:
            subplot = (4, 4, amp)
            ax = fig.add_subplot(*subplot)
            my_recarr = self.rec_array[np.where(self.rec_array.amp == amp)]
            plt.hist(my_recarr[pix0], color='blue', histtype='step',
                     range=(-10, 60), bins=20)
            plt.hist(my_recarr[pix1], color='red', histtype='step',
                     range=(-10, 60), bins=20)
            ax.set_xlabel('%(pix0)s, %(pix1)s (ADU)' % locals())
            ax.set_ylabel('entries / bin')
            ax.set_title('amp %i' % amp)
        return fig

    def pixel_diff_profile(self, pixel_coord='x', pix0='p3', pix1='p5',
                           figsize=(10, 10)):
        """
        Fit the difference of pix1 and pix0 profiles as a function of
        cluster peak pixel coordinate.
        """
        plt.rcParams['figure.figsize'] = figsize
        fig = plt.figure()
        data = []
        for amp in self.amps:
            subplot = (4, 4, amp)
            axes = fig.add_subplot(*subplot)
            my_recarr = self.rec_array[np.where(self.rec_array.amp == amp)]
            profile_plot(axes, my_recarr[pixel_coord], my_recarr[pix1],
                         color='red')
            profile_plot(axes, my_recarr[pixel_coord], my_recarr[pix0],
                         color='blue')
            prof = profile_plot(axes, my_recarr[pixel_coord],
                                (my_recarr[pix1] + my_recarr[pix0])/2.,
                                color='green')
            axes.set_title("amp %i" % amp)
            plt.xlim(min(self.rec_array[pixel_coord]),
                     max(self.rec_array[pixel_coord]))
            # Fit a line to the pix1 - pix0 profile.
            pars, cov = np.polyfit(prof['bin_centers'], prof['ymean'],
                                   1, cov=True)
            error = np.sqrt(cov[0][0])
            data.append([amp, pars[0], error, pars[0]/error])
        return fig, RecArray.create(data, names='amp slope error nsig'.split())

