"""
Module to manage plots for single sensor EO test reports.
"""
import os
import sys
import glob
import copy
import warnings
from collections import OrderedDict
import json
import numpy as np
import astropy.io.fits as fits
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.display.ds9 as ds9
import lsst.eotest.image_utils as imutils
from lsst.eotest.Estimator import Estimator
from lsst.eotest.sensor.ptcTask import ptc_func
from . import pylab_plotter as plot
from .MaskedCCD import MaskedCCD
from .EOTestResults import EOTestResults
from .Fe55GainFitter import Fe55GainFitter
from .fe55_psf import psf_sigma_statistics
from .DetectorResponse import DetectorResponse
from .crosstalk import CrosstalkMatrix
from .QE import QE_Data
from .AmplifierGeometry import parse_geom_kwd
from .cte_profile import *


class Subplot(object):
    def __init__(self, namps):
        self.nx = 4
        self.ny = namps/4

    def __call__(self, amp):
        return (self.nx, self.ny, amp)


def op_str(arg, format, op=None):
    if op is None:
        my_val = arg
    else:
        my_val = op(arg)
    try:
        return format % my_val
    except TypeError:
        return str(my_val)


def min_str(values, format):
    return op_str(values, format, op=np.min)


def max_str(values, format):
    return op_str(values, format, op=np.max)


def latex_minus_max(values, errors, format='%.2e'):
    # values or errors may contain nan's so eliminate those.
    index = np.where(~(np.isnan(values)) & ~(np.isnan(errors)))
    my_values = values[index]
    my_errors = errors[index]
    max_value = max(my_values)
    index = np.where(my_values == max_value)[0][0]
    if max_value < 0:
        template = '+ \\num{' + format + '}'
    else:
        template = '- \\num{' + format + '}'
    template += ' \pm \\num{' + format + '}'
    return template % (np.abs(max_value), errors[index])


def cmap_range(image_array, nsig=5):
    pixel_data = np.array(image_array, dtype=np.float).flatten()
    stats = afwMath.makeStatistics(pixel_data,
                                   afwMath.STDEVCLIP | afwMath.MEDIAN)
    median = stats.getValue(afwMath.MEDIAN)
    stdev = stats.getValue(afwMath.STDEVCLIP)
    vmin = max(min(pixel_data), median - nsig*stdev)
    vmax = min(max(pixel_data), median + nsig*stdev)
    return vmin, vmax

def flatten_across_segments(ccd, border=50):
    medians = dict()
    for amp in ccd:
        image = ccd.unbiased_and_trimmed_image(amp)
        bbox = image.getBBox().grow(-border)
        subimage = image.Factory(image, bbox)
        medians[amp] \
            = afwMath.makeStatistics(subimage, afwMath.MEDIAN).getValue()
    med0 = np.median(list(medians.values()))
    for amp in ccd:
        ccd[amp] *= med0/medians[amp]
    return ccd

def plot_flat(infile, nsig=3, cmap=pylab.cm.hot, win=None, subplot=(1, 1, 1),
              figsize=None, wl=None, gains=None, use_ds9=False, outfile=None,
              title=None, annotation='', flatten=False, binsize=1,
              bias_frame=None):
    ccd = MaskedCCD(infile, bias_frame=bias_frame)
    if flatten:
        ccd = flatten_across_segments(ccd)
    with fits.open(infile) as foo:
        if wl is None:
            # Extract wavelength from file
            try:
                wl = foo[0].header['MONOWL']
            except KeyError:
                wl = 0
        datasec = parse_geom_kwd(foo[1].header['DATASEC'])
        # Specialize to science sensor or wavefront sensor geometries.
        nx_segments = 8
        ny_segments = len(ccd)//nx_segments
        nx = nx_segments*(datasec['xmax'] - datasec['xmin'] + 1)
        ny = ny_segments*(datasec['ymax'] - datasec['ymin'] + 1)
        mosaic = np.zeros((ny, nx), dtype=np.float)
        amp_coords = {}
        for ypos in range(ny_segments):
            for xpos in range(nx_segments):
                amp = ypos*nx_segments + xpos + 1
                #
                # Determine subarray boundaries in the mosaicked image array
                # from DETSEC keywords for each segment.
                detsec = parse_geom_kwd(foo[amp].header['DETSEC'])
                xmin = nx - max(detsec['xmin'], detsec['xmax'])
                xmax = nx - min(detsec['xmin'], detsec['xmax']) + 1
                ymin = ny - max(detsec['ymin'], detsec['ymax'])
                ymax = ny - min(detsec['ymin'], detsec['ymax']) + 1
                if ny_segments == 1:
                    # We have a wavefront sensor, so only the "top"
                    # half of the sensor is available.
                    ymin += ny
                    ymax += ny
                #
                # Save coordinates of segment for later use
                #
                amp_coords[(xpos, ypos)] = xmin, xmax, ymin, ymax
                #
                # Extract the bias-subtracted masked image for this segment.
                segment_image = ccd.unbiased_and_trimmed_image(amp)
                subarr = segment_image.getImage().getArray()
                #
                # Determine flips in x- and y-directions in order to
                # get the (1, 1) pixel in the lower right corner.
                if detsec['xmax'] > detsec['xmin']:  # flip in x-direction
                    subarr = subarr[:, ::-1]
                if detsec['ymax'] > detsec['ymin']:  # flip in y-direction
                    subarr = subarr[::-1, :]
                #
                # Convert from ADU to e-
                if gains is not None:
                    subarr *= gains[amp]
                #
                # Set the subarray in the mosaicked image.
                mosaic[ymin:ymax, xmin:xmax] = subarr

    if binsize > 1:
        my_image = afwImage.ImageF(*mosaic.shape)
        my_image.array[:] = mosaic.transpose()
        mosaic = imutils.rebin(my_image, binsize=binsize).array/binsize**2
    #
    # Display the mosiacked image in ds9 using afwImage.
    if use_ds9:
        image = afwImage.ImageF(nx, ny)
        imarr = image.getArray()
        # This needs a flip in y to display properly in ds9 so that
        # amp 1 is in the lower right corner.
        imarr[:] = mosaic[::-1, :]
        ds9.mtv(image)
        ds9.ds9Cmd('zoom to fit')
        return
    #
    # Write a fits image with the mosaicked CCD data.
    if outfile is not None:
        hdulist = fits.HDUList()
        hdulist.append(fits.PrimaryHDU())
        hdulist[0].data = mosaic[::-1, :]
        warnings.filterwarnings('ignore', category=fits.verify.VerifyWarning,
                                append=True)
        hdulist.writeto(outfile, overwrite=True)
    #
    # Set the color map to extend over the range median +/- stdev(clipped)
    # of the pixel values.
    vmin, vmax = cmap_range(mosaic, nsig=nsig)
    if win is None:
        win = plot.Window(subplot=subplot, figsize=figsize,
                          xlabel='', ylabel='')
    else:
        win.select_subplot(*subplot)
    image = win.axes[-1].imshow(mosaic, interpolation='nearest', cmap=cmap)
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    image.set_norm(norm)
    if title is None:
        title = '%i nm' % wl
    win.axes[-1].set_title(title, fontsize='small')
    win.fig.colorbar(image)
    # Turn off tick labels for x- and y-axes
    pylab.setp(win.axes[-1].get_xticklabels(), visible=False)
    pylab.setp(win.axes[-1].get_yticklabels(), visible=False)
    # Label each segment by segment id.
    for ypos in range(ny_segments):
        for xpos in range(nx_segments):
            xmin, xmax, ymin, ymax = amp_coords[(xpos, ypos)]
            xx = float(xmax + xmin)/2./float(nx)
            yy = 1. - (float(ymax)/float(ny) - 0.05)
            if yy > 0.5:
                yy = 1 - (yy - 0.5)
            amp = ypos*nx_segments + xpos + 1
            if ny_segments == 1:
                # We have a wavefront sensor which just has the
                # "lower" half of the CCD, so shift the amp number by
                # 8 segments.
                amp += 8
            seg_id = imutils.channelIds[amp]
            pylab.annotate('%s' % seg_id, (xx, yy), xycoords='axes fraction',
                           size='x-small', horizontalalignment='center')
    plt.annotate(annotation, (1, -0.1), xycoords='axes fraction',
                 horizontalalignment='right', verticalalignment='bottom')
    return win


def fe55_zoom(infile, size=250, amp=1, cmap=pylab.cm.hot, nsig=10,
              subplot=(1, 1, 1), win=None, figsize=None, title=None,
              axisRange=None, annotation=''):
    ccd = MaskedCCD(infile)
    image = ccd[amp].getImage()
    nymax, nxmax = image.getArray().shape
    sub_image = afwImage.ImageF(size, size)
    sub_image.getArray()[:] = image.getArray()[:-size-1:-1, -size:]
    if win is None:
        win = plot.Window(subplot=subplot, figsize=figsize,
                          xlabel='', ylabel='')
    else:
        win.select_subplot(*subplot)
    if axisRange is None:
        axisRange = (nxmax-size, nxmax, nymax-size, nymax)
    image = win.axes[-1].imshow(sub_image.getArray(), interpolation='nearest',
                                cmap=cmap, extent=axisRange)
    pylab.xlabel('x pixel')
    pylab.ylabel('y pixel')
    vmin, vmax = cmap_range(sub_image.getArray(), nsig=nsig)
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    image.set_norm(norm)
    if title is None:
        title = os.path.basename(infile)
    win.axes[-1].set_title(title)
    win.fig.colorbar(image)
    plt.annotate(annotation, (1, -0.1), xycoords='axes fraction',
                 horizontalalignment='right', verticalalignment='bottom')
    return win

class OverscanTestPlots(object):

    def __init__(self, sensor_id, overscan_file=None, title_addendum=None):

        self.sensor_id = sensor_id
        self._title_addendum = title_addendum
        self._plot_title = None
        if overscan_file is None:
            overscan_file = self._fullpath('%s_overscan_results.fits' % sensor_id)
        if not os.path.exists(overscan_file):
            raise RuntimeError("OverscanTestPlots: %s not found" % overscan_file)

        with fits.open(overscan_file) as hdul:
            self.header = hdul[0].header
            self.all_amps = range(1, int(self.header['NAMPS'])+1)
            self.overscan_data = {i : hdul[i].data for i in self.all_amps}

    @property
    def plot_title(self):
        if self._plot_title is None:
            self._plot_title = self.sensor_id
            if self._title_addendum is not None:
                self._plot_title \
                    = ' '.join((self._plot_title, self._title_addendum))
        return self._plot_title

    def _fullpath(self, basename):
        return os.path.join(self.rootdir, basename)

    def parallel_cti_curves(self, figsize=(12, 8)):

        fig, ax = plt.subplots(1, 1, figsize=figsize)

        datasec = self.header['DATASEC']
        amp_geom = parse_geom_kwd(datasec)
        num_transfers = amp_geom['ymax']

        cmap = plt.get_cmap("tab10")

        for amp in self.all_amps:

            if amp > 8: marker='s'
            else: marker = '^'

            data = self.overscan_data[amp]['ROW_MEAN']
            signal = self.overscan_data[amp]['FLATFIELD_SIGNAL']
            lastpixel = data[:, 0]
            overscan1 = data[:, 1]
            overscan2 = data[:, 2]
            cti = (overscan1+overscan2)/(num_transfers*lastpixel)

            indices = signal<175000.

            ax.plot(signal[indices], cti[indices], label='Amp {0}'.format(amp),
                    color = cmap((amp-1)%8), marker=marker, markersize=5)

        ax.axhline(y=3.E-6, linestyle='--', color='black')

        ax.tick_params(axis='x', labelsize=14)
        ax.tick_params(axis='y', labelsize=14)
        ax.tick_params(axis='both', which='major', length=8, width=1)
        ax.tick_params(axis='both', which='minor', length=4, width=1)
        ax.set_ylim(bottom=5E-8, top=2E-4)
        ax.set_xlim(left=50.0, right=240000.)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.grid(True, which='major', axis='both')
        ax.set_xlabel('Flat Field Signal [e-]', fontsize=18)
        ax.set_ylabel('Parallel CTI', fontsize=18)
        ax.legend(fontsize=12, loc=1, ncol=4)
        ax.set_title('Parallel CTI from EPER, {0}'.format(self.plot_title),
                     fontsize=18)

    def parallel_eper_high_curves(self, figsize=(12, 8)):

        fig, axes = plt.subplots(4, 4, sharey=True, sharex=True, figsize=figsize)
        axes = axes.flatten()

        target_signals = [25000, 50000, 75000, 100000, 150000]

        for amp in self.all_amps:

            signal = self.overscan_data[amp]['FLATFIELD_SIGNAL']

            for target in target_signals:

                i = min(range(signal.shape[0]),
                        key=lambda i: abs(signal[i]-target))

                data = self.overscan_data[amp]['ROW_MEAN'][i, :]
                oscan = data[1:]
                columns = np.arange(1, data.shape[0])

                axes[amp-1].plot(columns, oscan,
                                 label='{0:.0f} ke-'.format(target/1000.))

            axes[amp-1].set_yscale('log')
            axes[amp-1].set_ylim(.05, 400)
            axes[amp-1].tick_params(axis='x', labelsize=12)
            axes[amp-1].tick_params(axis='y', labelsize=12)
            axes[amp-1].set_xlim(0., 30.)
            axes[amp-1].grid(True, which='major', axis='both')
            axes[amp-1].set_title('Amp {0}'.format(amp), fontsize=12)
            axes[amp-1].tick_params(axis='both', which='minor')

        h, l = axes[-1].get_legend_handles_labels()
        fig.subplots_adjust(bottom=0.12)
        fig.legend(h, l, loc='lower center', ncol=len(target_signals),
                   fontsize=14)

        fig.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor='none', top=False, bottom=False,
                        left=False, right=False)
        plt.xlabel('Overscan Row Number', fontsize=14)
        plt.ylabel('Mean Row Signal [e-]', fontsize=14, labelpad=15)
        plt.title('High Signal Parallel EPER, {0}'.format(self.plot_title),
                  fontsize=16, pad=20)

    def parallel_eper_low_curves(self, figsize=(12, 8)):

        fig, axes = plt.subplots(4, 4, sharey=True, sharex=True, figsize=figsize)
        axes = axes.flatten()

        target_signals = [100, 500, 1000, 2500, 5000]
        for amp in self.all_amps:

            signal = self.overscan_data[amp]['FLATFIELD_SIGNAL']

            for target in target_signals:

                i = min(range(signal.shape[0]),
                        key=lambda i: abs(signal[i]-target))

                data = self.overscan_data[amp]['ROW_MEAN'][i, :]
                oscan = data[1:]
                columns = np.arange(1, data.shape[0])

                axes[amp-1].plot(columns, oscan,
                                 label='{0:d} e-'.format(target))

            axes[amp-1].set_xlim(0.5, 5.5)
            axes[amp-1].tick_params(axis='x', labelsize=12)
            axes[amp-1].tick_params(axis='y', labelsize=12)
            axes[amp-1].grid(True, which='major', axis='both')
            axes[amp-1].set_title('Amp {0}'.format(amp), fontsize=14)


            if (amp-1) % 4 == 0:
                axes[amp-1].set_ylim(-4, 22)
                axes[amp-1].set_yticks([0, 5, 10, 15, 20])

        h, l = axes[-1].get_legend_handles_labels()
        fig.subplots_adjust(bottom=0.12)
        fig.legend(h, l, loc='lower center', ncol=len(target_signals),
                   fontsize=14)

        fig.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor='none', top=False, bottom=False,
                        left=False, right=False)
        plt.xlabel('Overscan Row Number', fontsize=14)
        plt.ylabel('Mean Row Signal [e-]', fontsize=14)
        plt.title('Low Signal Parallel EPER, {0}'.format(self.plot_title),
                  fontsize=16, pad=20)

    def parallel_overscan_noise_curves(self, figsize=(12, 8)):

        fig, ax = plt.subplots(1, 1, figsize=figsize)

        cmap = plt.get_cmap("tab10")

        for amp in self.all_amps:

            if amp > 8: marker='s'
            else: marker = '^'

            signal = self.overscan_data[amp]['FLATFIELD_SIGNAL']
            noise = self.overscan_data[amp]['PARALLEL_OVERSCAN_NOISE']

            ax.plot(signal, noise, label='Amp {0}'.format(amp), linestyle='-',
                    marker=marker, color = cmap((amp-1)%8), markersize=5)

        ax.set_xscale('log')
        ax.set_xlim(left=50, right=240000)
        ax.set_ylim(bottom=0.0, top=min(21.0, np.max(noise)+2.0))
        ax.grid(True, which='major', axis='both')
        ax.set_xlabel('Flat Field Signal [e-]', fontsize=18)
        ax.set_ylabel('Parallel Overscan Noise [e-]', fontsize=18)
        ax.tick_params(axis='x', labelsize=14)
        ax.tick_params(axis='y', labelsize=14)

        plt.legend(fontsize=14,  loc = 'upper left', ncol=4)
        plt.title('Parallel Overscan Pixel Noise, {0}'.format(self.plot_title),
                  fontsize=18)

    def parallel_overscan_signal_curves(self, figsize=(12, 8)):

        fig, axes = plt.subplots(4, 4, sharey=True, sharex=True, figsize=figsize)
        axes = axes.flatten()

        for amp in self.all_amps:

            data = self.overscan_data[amp]['ROW_MEAN']
            signal = self.overscan_data[amp]['FLATFIELD_SIGNAL']
            overscan1 = data[:, 1]
            overscan2 = data[:, 2]

            axes[amp-1].plot(signal, overscan1, label='Overscan 1')
            axes[amp-1].plot(signal, overscan2, label='Overscan 2')

            axes[amp-1].set_yscale('symlog', threshold=1.0)
            axes[amp-1].set_xscale('log')
            axes[amp-1].set_ylim(bottom=-1.0, top=2000.)
            axes[amp-1].set_xlim(left=50, right=240000)
            axes[amp-1].set_title('Amp {0}'.format(amp), fontsize=12)
            axes[amp-1].grid(True, which='major', axis='both')

        h, l = axes[-1].get_legend_handles_labels()
        fig.subplots_adjust(bottom=0.12)
        fig.legend(h, l, loc='lower center', ncol=2, fontsize=14)

        fig.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor='none', top=False, bottom=False,
                        left=False, right=False)
        plt.ylabel('Parallel Overscan Mean Signal [e-]', fontsize=14, labelpad=15)
        plt.xlabel('Flat Field Signal [e-]', fontsize=14)
        plt.title('Overscan Pixel Signal, {0}'.format(self.plot_title),
                  fontsize=16, pad=20)

    def parallel_overscan_sum_curves(self, figsize=(12, 8)):

        fig, ax = plt.subplots(1, 1, figsize=figsize)

        cmap = plt.get_cmap("tab10")

        for amp in self.all_amps:

            if amp > 8: marker = 's'
            else: marker = '^'

            data = self.overscan_data[amp]['ROW_MEAN']
            oscansum = np.sum(data[:, 3:20], axis=1)
            signal = self.overscan_data[amp]['FLATFIELD_SIGNAL']

            ax.plot(signal, oscansum, label='Amp {0}'.format(amp),
                    linestyle='-', marker=marker, color = cmap((amp-1)%8),
                    markersize=5)

        ax.set_xscale('log')
        ax.set_xlim(left=50.0, right=240000.0)
        ax.set_ylim(bottom=-7.0, top=52.0)
        ax.grid(True, which='major', axis='both')
        ax.set_xlabel('Flat Field Signal [e-]', fontsize=18)
        ax.set_ylabel('Parallel Overscan Mean Sum [e-]', fontsize=18)
        ax.tick_params(axis='x', labelsize=14)
        ax.tick_params(axis='y', labelsize=14)

        ax.legend(fontsize=14,  loc = 'upper left', ncol=4)
        ax.set_title('Summed Parallel Overscan Pixel Signal [3:20], {0}'.format(self.plot_title),
                     fontsize=18)

    def serial_cti_curves(self, figsize=(12, 8)):

        fig, ax = plt.subplots(1, 1, figsize=figsize)

        datasec = self.header['DATASEC']
        amp_geom = parse_geom_kwd(datasec)
        num_transfers = amp_geom['xmax']

        cmap = plt.get_cmap("tab10")

        for amp in self.all_amps:

            if amp > 8: marker='s'
            else: marker = '^'

            data = self.overscan_data[amp]['COLUMN_MEAN']
            signal = self.overscan_data[amp]['FLATFIELD_SIGNAL']
            lastpixel = data[:, 0]
            overscan1 = data[:, 1]
            overscan2 = data[:, 2]
            cti = (overscan1+overscan2)/(num_transfers*lastpixel)

            indices = signal<175000.

            ax.plot(signal[indices], cti[indices], label='Amp {0}'.format(amp),
                    color = cmap((amp-1)%8), marker=marker, markersize=5)

        ax.axhline(y=5.E-6, linestyle='--', color='black')

        ax.tick_params(axis='x', labelsize=14)
        ax.tick_params(axis='y', labelsize=14)
        ax.tick_params(axis='both', which='major', length=8, width=1)
        ax.tick_params(axis='both', which='minor', length=4, width=1)
        ax.set_ylim(bottom=7E-8, top=2E-4)
        ax.set_xlim(left=50.0, right=240000.)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.grid(True, which='major', axis='both')
        ax.set_xlabel('Flat Field Signal [e-]', fontsize=18)
        ax.set_ylabel('Serial CTI', fontsize=18)
        ax.legend(fontsize=12, loc=1, ncol=4)
        ax.set_title('Serial CTI from EPER, {0}'.format(self.plot_title),
                     fontsize=18)

    def serial_eper_high_curves(self, figsize=(12, 8)):

        fig, axes = plt.subplots(4, 4, sharey=True, sharex=True, figsize=figsize)
        axes = axes.flatten()

        target_signals = [25000, 50000, 75000, 100000, 150000]

        for amp in self.all_amps:

            signal = self.overscan_data[amp]['FLATFIELD_SIGNAL']

            for target in target_signals:

                i = min(range(signal.shape[0]),
                        key=lambda i: abs(signal[i]-target))

                data = self.overscan_data[amp]['COLUMN_MEAN'][i, :]
                oscan = data[1:]
                columns = np.arange(1, data.shape[0])

                axes[amp-1].plot(columns, oscan,
                                 label='{0:.0f} ke-'.format(target/1000.))

            axes[amp-1].set_yscale('log')
            axes[amp-1].set_ylim(.05, 400)
            axes[amp-1].tick_params(axis='x', labelsize=12)
            axes[amp-1].tick_params(axis='y', labelsize=12)
            axes[amp-1].set_xlim(0., 30.)
            axes[amp-1].grid(True, which='major', axis='both')
            axes[amp-1].set_title('Amp {0}'.format(amp), fontsize=12)
            axes[amp-1].tick_params(axis='both', which='minor')

        h, l = axes[-1].get_legend_handles_labels()
        fig.subplots_adjust(bottom=0.12)
        fig.legend(h, l, loc='lower center', ncol=len(target_signals),
                   fontsize=14)

        fig.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor='none', top=False, bottom=False,
                        left=False, right=False)
        plt.xlabel('Overscan Column Number', fontsize=14)
        plt.ylabel('Mean Column Signal [e-]', fontsize=14, labelpad=15)
        plt.title('High Signal Serial EPER, {0}'.format(self.plot_title),
                  fontsize=16, pad=20)

    def serial_eper_low_curves(self, figsize=(12, 8)):

        fig, axes = plt.subplots(4, 4, sharey=True, sharex=True, figsize=figsize)
        axes = axes.flatten()

        target_signals = [100, 500, 1000, 2500, 5000]
        for amp in self.all_amps:

            signal = self.overscan_data[amp]['FLATFIELD_SIGNAL']

            for target in target_signals:

                i = min(range(signal.shape[0]),
                        key=lambda i: abs(signal[i]-target))

                data = self.overscan_data[amp]['COLUMN_MEAN'][i, :]
                oscan = data[1:]
                columns = np.arange(1, data.shape[0])

                axes[amp-1].plot(columns, oscan,
                                 label='{0:d} e-'.format(target))

            axes[amp-1].set_xlim(0.5, 5.5)
            axes[amp-1].tick_params(axis='x', labelsize=12)
            axes[amp-1].tick_params(axis='y', labelsize=12)
            axes[amp-1].grid(True, which='major', axis='both')
            axes[amp-1].set_title('Amp {0}'.format(amp), fontsize=14)


            if (amp-1) % 4 == 0:
                axes[amp-1].set_ylim(-4, 22)
                axes[amp-1].set_yticks([0, 5, 10, 15, 20])

        h, l = axes[-1].get_legend_handles_labels()
        fig.subplots_adjust(bottom=0.12)
        fig.legend(h, l, loc='lower center', ncol=len(target_signals),
                   fontsize=14)

        fig.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor='none', top=False, bottom=False,
                        left=False, right=False)
        plt.xlabel('Overscan Column Number', fontsize=14)
        plt.ylabel('Mean Column Signal [e-]', fontsize=14)
        plt.title('Low Signal Serial EPER, {0}'.format(self.plot_title),
                  fontsize=16, pad=20)

    def serial_overscan_noise_curves(self, figsize=(12, 8)):

        fig, ax = plt.subplots(1, 1, figsize=figsize)

        cmap = plt.get_cmap("tab10")

        for amp in self.all_amps:

            if amp > 8: marker='s'
            else: marker = '^'

            signal = self.overscan_data[amp]['FLATFIELD_SIGNAL']
            noise = self.overscan_data[amp]['SERIAL_OVERSCAN_NOISE']

            ax.plot(signal, noise, label='Amp {0}'.format(amp), linestyle='-',
                    marker=marker, color = cmap((amp-1)%8), markersize=5)

        ax.set_xscale('log')
        ax.set_xlim(left=50, right=240000)
        ax.set_ylim(bottom=0.0, top=min(21.0, np.max(noise)+2.0))
        ax.grid(True, which='major', axis='both')
        ax.set_xlabel('Flat Field Signal [e-]', fontsize=18)
        ax.set_ylabel('Serial Overscan Noise [e-]', fontsize=18)
        ax.tick_params(axis='x', labelsize=14)
        ax.tick_params(axis='y', labelsize=14)

        plt.legend(fontsize=14,  loc = 'upper left', ncol=4)
        plt.title('Serial Overscan Pixel Noise, {0}'.format(self.plot_title),
                  fontsize=18)

    def serial_overscan_signal_curves(self, figsize=(12, 8)):

        fig, axes = plt.subplots(4, 4, sharey=True, sharex=True, figsize=figsize)
        axes = axes.flatten()

        for amp in self.all_amps:

            data = self.overscan_data[amp]['COLUMN_MEAN']
            signal = self.overscan_data[amp]['FLATFIELD_SIGNAL']
            overscan1 = data[:, 1]
            overscan2 = data[:, 2]

            axes[amp-1].plot(signal, overscan1, label='Overscan 1')
            axes[amp-1].plot(signal, overscan2, label='Overscan 2')

            axes[amp-1].set_yscale('symlog', threshold=1.0)
            axes[amp-1].set_xscale('log')
            axes[amp-1].set_ylim(bottom=-1.0, top=2000.)
            axes[amp-1].set_xlim(left=50, right=240000)
            axes[amp-1].set_title('Amp {0}'.format(amp), fontsize=12)
            axes[amp-1].grid(True, which='major', axis='both')

        h, l = axes[-1].get_legend_handles_labels()
        fig.subplots_adjust(bottom=0.12)
        fig.legend(h, l, loc='lower center', ncol=2, fontsize=14)

        fig.add_subplot(111, frameon=False)
        plt.tick_params(labelcolor='none', top=False, bottom=False,
                        left=False, right=False)
        plt.ylabel('Serial Overscan Mean Signal [e-]', fontsize=14, labelpad=15)
        plt.xlabel('Flat Field Signal [e-]', fontsize=14)
        plt.title('Overscan Pixel Signal, {0}'.format(self.plot_title),
                  fontsize=16, pad=20)

    def serial_overscan_sum_curves(self, figsize=(12, 8)):

        fig, ax = plt.subplots(1, 1, figsize=figsize)

        cmap = plt.get_cmap("tab10")

        for amp in self.all_amps:

            if amp > 8: marker = 's'
            else: marker = '^'

            data = self.overscan_data[amp]['COLUMN_MEAN']
            oscansum = np.sum(data[:, 5:25], axis=1)
            signal = self.overscan_data[amp]['FLATFIELD_SIGNAL']

            ax.plot(signal, oscansum, label='Amp {0}'.format(amp),
                    linestyle='-', marker=marker, color = cmap((amp-1)%8),
                    markersize=5)

        ax.set_xscale('log')
        ax.set_xlim(left=50.0, right=240000.0)
        ax.set_ylim(bottom=-7.0, top=52.0)
        ax.grid(True, which='major', axis='both')
        ax.set_xlabel('Flat Field Signal [e-]', fontsize=18)
        ax.set_ylabel('Serial Overscan Mean Sum [e-]', fontsize=18)
        ax.tick_params(axis='x', labelsize=14)
        ax.tick_params(axis='y', labelsize=14)

        ax.legend(fontsize=14,  loc = 'upper left', ncol=4)
        ax.set_title('Summed Serial Overscan Pixel Signal [5:25], {0}'.format(self.plot_title),
                     fontsize=18)

class EOTestPlots(object):
    band_pass = QE_Data.band_pass
    prnu_wls = (350, 450, 500, 620, 750, 870, 1000)

    def __init__(self, sensor_id, rootdir='.', output_dir='.',
                 interactive=False, results_file=None, xtalk_file=None,
                 title_addendum=None):
        self.sensor_id = sensor_id
        self.rootdir = rootdir
        self.output_dir = output_dir
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        self.interactive = interactive
        plot.pylab.interactive(interactive)
        self._title_addendum = title_addendum
        self._plot_title = None
        if results_file is None:
            results_file = self._fullpath('%s_eotest_results.fits' % sensor_id)
        if not os.path.exists(results_file):
            raise RuntimeError("EOTestPlots: %s not found" % results_file)
        self.results = EOTestResults(results_file)
        self._qe_data = None
        self._qe_file = self._fullpath('%s_QE.fits' % self.sensor_id)
        self.specs = CcdSpecs(results_file, plotter=self,
                              xtalk_file=xtalk_file, prnu_wls=self.prnu_wls)
        self._linearity_results = None
        self.subplot = Subplot(len(self.results['AMP']))

    @property
    def plot_title(self):
        if self._plot_title is None:
            self._plot_title = self.sensor_id
            if self._title_addendum is not None:
                self._plot_title \
                    = ' '.join((self._plot_title, self._title_addendum))
        return self._plot_title

    @property
    def qe_data(self):
        if self._qe_data is None:
            self._qe_data = fits.open(self._qe_file)
        return self._qe_data

    def _save_fig(self, outfile_root):
        plot.pylab.savefig(self._outputpath('%s.png' % outfile_root))

    def _fullpath(self, basename):
        return os.path.join(self.rootdir, basename)

    def _outputpath(self, basename):
        return os.path.join(self.output_dir, basename)

    def crosstalk_matrix(self, cmap=pylab.cm.hot, xtalk_file=None):
        if xtalk_file is None:
            xtalk_file = os.path.join(self.rootdir,
                                      '%s_xtalk_matrix.fits' % self.sensor_id)
        foo = CrosstalkMatrix(xtalk_file)
        foo.plot(title="Crosstalk, %s" % self.plot_title)
        return foo

    def persistence(self, infile=None, figsize=(11, 8.5)):
        if infile is None:
            infile = self._fullpath('%s_persistence.fits' % self.sensor_id)
        results = fits.open(infile)
        times = results[1].data.field('TIME')
        win = None
        all_amps = imutils.allAmps(infile)
        for amp in all_amps:
            flux = results[1].data.field('MEDIAN%02i' % amp)
            stdev = results[1].data.field('STDEV%02i' % amp)
            subplot = self.subplot(amp)
            if amp == 1:
                win = plot.Window(subplot=subplot, figsize=figsize,
                                  xlabel=r'Time since end of flat exposure (s)',
                                  ylabel=r'Deferred charge (e-/pixel)')
                win.frameAxes.text(0.5, 1.08,
                                   'Image Persistence vs Time, %s'
                                   % self.sensor_id,
                                   horizontalalignment='center',
                                   verticalalignment='top',
                                   transform=win.frameAxes.transAxes,
                                   size='large')
            else:
                win.select_subplot(*subplot)
            self._offset_subplot(win)
            try:
                plot.xyplot(times, flux, yerr=stdev, xname='', yname='',
                            new_win=False)
                pylab.annotate('Amp %i' % amp, (0.5, 0.8),
                               xycoords='axes fraction', size='x-small')
            except Exception as eobj:
                print("Exception raised in generating image persistence plot for amp", amp)
                print(eobj)
                # Continue with remaining amps

    def psf_dists(self, chiprob_min=0.1, fe55_file=None, figsize=(11, 8.5),
                  xrange=(2, 6), bins=50):
        if fe55_file is None:
            fe55_file = glob.glob(self._fullpath('%s_psf_results*.fits'
                                                 % self.sensor_id))[0]
        fe55_catalog = fits.open(fe55_file)
        win = None
        for amp in imutils.allAmps(fe55_file):
            subplot = self.subplot(amp)
            chiprob = fe55_catalog[amp].data.field('CHIPROB')
            index = np.where(chiprob > chiprob_min)
            # Read sigma values and convert from pixels to microns
            sigmax = fe55_catalog[amp].data.field('SIGMAX')[index]*10.
            sigmay = fe55_catalog[amp].data.field('SIGMAY')[index]*10.
            sigma = sorted(np.concatenate((sigmax, sigmay)))
            if amp == 1:
                win = plot.Window(subplot=subplot, figsize=figsize,
                                  xlabel=r'PSF sigma ($\mu$); $\sigma_x$ in blue, $\sigma_y$ in red',
                                  ylabel=r'entries / bin', size='large')
                win.frameAxes.text(0.5, 1.08,
                                   'PSF from Fe55 data, %s' % self.sensor_id,
                                   horizontalalignment='center',
                                   verticalalignment='top',
                                   transform=win.frameAxes.transAxes,
                                   size='large')
            else:
                win.select_subplot(*subplot)
            self._offset_subplot(win)
            try:
                plot.histogram(sigmax, xrange=xrange, bins=bins, new_win=False,
                               xname='', yname='', color='blue')
                plot.histogram(sigmay, oplot=1, color='red', xrange=xrange,
                               bins=bins)
                plot.histogram(sigma, oplot=1, xrange=xrange, bins=bins)
                pylab.xticks(list(range(xrange[0], xrange[1]+1)))
                # Find mode from histogram data
                mode, median, mean = psf_sigma_statistics(sigma, bins=bins,
                                                          range=xrange)
                plot.vline(5)
                plot.vline(mode, color='r')
                pylab.annotate('Amp %i\nmode=%.2f' % (amp, mode), (0.5, 0.8),
                               xycoords='axes fraction', size='x-small')
            except Exception as eobj:
                # Skip this plot so that the rest of the plots can be
                # generated.
                print("Exception raised in generating PSF sigma plot for amp", amp)
                print(eobj)

    def fe55_dists(self, chiprob_min=0.1, fe55_file=None, figsize=(11, 8.5),
                   hist_nsig=10, xrange_scale=1, dn_range=None):
        if fe55_file is None:
            fe55_file = glob.glob(self._fullpath('%s_psf_results*.fits'
                                                 % self.sensor_id))[0]
        fe55_catalog = fits.open(fe55_file)
        win = None
        for amp in imutils.allAmps(fe55_file):
            chiprob = fe55_catalog[amp].data.field('CHIPROB')
            index = np.where(chiprob > chiprob_min)
            dn = fe55_catalog[amp].data.field('DN')[index]
            if dn_range is not None:
                dn = dn[np.where((dn >= dn_range[0]) & (dn <= dn_range[1]))]
            foo = Fe55GainFitter(dn)
            try:
                foo.fit(hist_nsig=hist_nsig)
            except:
                continue
            if win is None:
                win = foo.plot(interactive=self.interactive,
                               subplot=self.subplot(amp),
                               figsize=figsize, frameLabels=True, amp=amp,
                               xrange_scale=xrange_scale)
                win.frameAxes.text(0.5, 1.08, 'Fe55, %s' % self.sensor_id,
                                   horizontalalignment='center',
                                   verticalalignment='top',
                                   transform=win.frameAxes.transAxes,
                                   size='large')
            else:
                foo.plot(interactive=self.interactive,
                         subplot=self.subplot(amp), win=win,
                         frameLabels=True, amp=amp, xrange_scale=xrange_scale)
            pylab.locator_params(axis='x', nbins=4, tight=True)

    def ptcs(self, xrange=None, yrange=None, figsize=(11, 8.5), ptc_file=None):
        if ptc_file is None:
            ptc_file = self._fullpath('%s_ptc.fits' % self.sensor_id)
        with fits.open(ptc_file) as ptc:
            for amp in imutils.allAmps(ptc_file):
                mean = ptc[1].data.field('AMP%02i_MEAN' % amp)
                var = ptc[1].data.field('AMP%02i_VAR' % amp)
                subplot = self.subplot(amp)
                if amp == 1:
                    win = plot.Window(subplot=subplot, figsize=figsize,
                                      xlabel=r'mean (ADU)',
                                      ylabel=r'variance (ADU$^2$)', size='large')
                    win.frameAxes.text(0.5, 1.08,
                                       'Photon Transfer Curves, %s'
                                       % self.sensor_id,
                                       horizontalalignment='center',
                                       verticalalignment='top',
                                       transform=win.frameAxes.transAxes,
                                       size='large')
                else:
                    win.select_subplot(*subplot)
                self._offset_subplot(win)
                win = plot.xyplot(mean, var, xname='', yname='',
                                  xrange=xrange, yrange=yrange,
                                  xlog=1, ylog=1, new_win=False,)
                axes = pylab.gca()
                xx = np.logspace(np.log10(min(mean)), np.log10(max(mean)), 20)
                # Plot PTC curves using gain measurements.
                ptc_gain = self.results['PTC_GAIN'][amp-1]
                ptc_gain_error = self.results['PTC_GAIN_ERROR'][amp-1]
                ptc_noise = self.results['PTC_NOISE'][amp-1]
                ptc_a00 = self.results['PTC_A00'][amp-1]
                ptc_a00_error = self.results['PTC_A00_ERROR'][amp-1]
                ptc_turnoff = self.results['PTC_TURNOFF'][amp-1]
                plot.curve(xx, ptc_func((ptc_a00, ptc_gain, ptc_noise*ptc_noise), xx),
                           oplot=1, color='b', lineStyle=':')
                note = 'Amp %i\nGain = %.2f +/- %.2f\nA00 = %.1e +/- %.1e \nTurnoff = %7.0f'\
                    % (amp, ptc_gain, ptc_gain_error, ptc_a00, ptc_a00_error,
                       ptc_turnoff)
                pylab.annotate(note, (0.05, 0.9), xycoords='axes fraction',
                               verticalalignment='top', size='x-small')

    def bf_curves(self, xrange=None, yrange=None, figsize=(24, 12),
                  bf_file=None, adu_max=1e5):
        if bf_file is None:
            bf_file = self._fullpath('%s_bf.fits' % self.sensor_id)
        fig = plt.figure(figsize=figsize)
        with fits.open(bf_file) as bf:

            fig.add_subplot(3, 2, 1)
            for amp in imutils.allAmps(bf_file):
                mean = bf[1].data.field('AMP%02i_MEAN' % amp)
                xcorr = bf[1].data.field('AMP%02i_COV10' % amp)
                index = np.argsort(mean)
                mean = mean[index]
                xcorr = xcorr[index]
                index = np.where(mean < adu_max)
                plt.plot(mean[index], xcorr[index]/mean[index], label='%s' % amp)
            plt.xlabel('mean signal (ADU)', fontsize='small')
            plt.ylabel('cov(1, 0)/mean', fontsize='small')
            plt.legend(fontsize='x-small', loc=2)
            plt.title('Brighter-Fatter cov(1, 0), %s' % self.plot_title,
                      fontsize='small')

            fig.add_subplot(3, 2, 2)
            for amp in imutils.allAmps(bf_file):
                mean = bf[1].data.field('AMP%02i_MEAN' % amp)
                xcorr = bf[1].data.field('AMP%02i_COV20' % amp)
                index = np.argsort(mean)
                mean = mean[index]
                xcorr = xcorr[index]
                index = np.where(mean < adu_max)
                plt.plot(mean[index], xcorr[index]/mean[index], label='%s' % amp)
            plt.xlabel('mean signal (ADU)', fontsize='small')
            plt.ylabel('cov(2, 0)/mean', fontsize='small')
            plt.legend(fontsize='x-small', loc=2)
            plt.title('Brighter-Fatter cov(2, 0), %s' % self.plot_title,
                      fontsize='small')

            fig.add_subplot(3, 2, 3)
            for amp in imutils.allAmps(bf_file):
                mean = bf[1].data.field('AMP%02i_MEAN' % amp)
                ycorr = bf[1].data.field('AMP%02i_COV01' % amp)
                index = np.argsort(mean)
                mean = mean[index]
                ycorr = ycorr[index]
                index = np.where(mean < adu_max)
                plt.plot(mean[index], ycorr[index]/mean[index],
                         label='%s' % amp)
            plt.xlabel('mean signal (ADU)', fontsize='small')
            plt.ylabel('cov(0, 1)/mean', fontsize='small')
            plt.legend(fontsize='x-small', loc=2)
            plt.title('Brighter-Fatter cov(0, 1), %s' % self.plot_title,
                      fontsize='small')
            fig.add_subplot(3, 2, 4)
            for amp in imutils.allAmps(bf_file):
                mean = bf[1].data.field('AMP%02i_MEAN' % amp)
                ycorr = bf[1].data.field('AMP%02i_COV02' % amp)
                index = np.argsort(mean)
                mean = mean[index]
                ycorr = ycorr[index]
                index = np.where(mean < adu_max)
                plt.plot(mean[index], ycorr[index]/mean[index], label='%s' % amp)
            plt.xlabel('mean signal (ADU)', fontsize='small')
            plt.ylabel('cov(0, 2)/mean', fontsize='small')
            plt.legend(fontsize='x-small', loc=2)
            plt.title('Brighter-Fatter cov(0, 2), %s' % self.plot_title,
                      fontsize='small')

            fig.add_subplot(3, 2, 5)
            for amp in imutils.allAmps(bf_file):
                mean = bf[1].data.field('AMP%02i_MEAN' % amp)
                ycorr = bf[1].data.field('AMP%02i_COV11' % amp)
                index = np.argsort(mean)
                mean = mean[index]
                ycorr = ycorr[index]
                index = np.where(mean < adu_max)
                plt.plot(mean[index], ycorr[index]/mean[index], label='%s' % amp)
            plt.xlabel('mean signal (ADU)', fontsize='small')
            plt.ylabel('cov(1, 1)/mean', fontsize='small')
            plt.legend(fontsize='x-small', loc=2)
            plt.title('Brighter-Fatter cov(1, 1), %s' % self.plot_title,
                      fontsize='small')
        plt.tight_layout()

    def _offset_subplot(self, win, xoffset=0.025, yoffset=0.025):
        bbox = win.axes[-1].get_position()
        points = bbox.get_points()
        points[0] += xoffset
        points[1] += yoffset
        bbox.set_points(points)
        win.axes[-1].set_position(bbox)

    def gains(self, oplot=0, xoffset=0.25, width=0.5, xrange=None):
        results = self.results
        gain = results['GAIN']
        error = results['GAIN_ERROR']
        try:
            ptc_gain = results['PTC_GAIN']
            ptc_error = results['PTC_GAIN_ERROR']
            ymin = min(max(min(gain - error), min(gain - 1)),
                       max(min(ptc_gain - ptc_error), min(ptc_gain - 1)))
            ymax = max(min(max(gain + error), max(gain + 1)),
                       min(max(ptc_gain + ptc_error), max(ptc_gain + 1)))
            yname = 'gain (e-/DN) (Fe55: black; PTC: red)'
        except KeyError:
            ymin = max(min(gain - error), min(gain - 1))
            ymax = min(max(gain + error), max(gain + 1))
            yname = 'gain (e-/DN)'
        if xrange is None:
            xrange = (0, len(gain) + 0.5)
        win = plot.xyplot(results['AMP'], results['GAIN'],
                          yerr=results['GAIN_ERROR'], xname='AMP',
                          yname=yname, xrange=xrange, yrange=(ymin, ymax))
        try:
            plot.xyplot(results['AMP'], results['PTC_GAIN'],
                        yerr=results['PTC_GAIN_ERROR'], oplot=1, color='r')
        except:
            pass
        win.set_title("System Gain, %s" % self.plot_title)

    def ptc_bf(self, oplot=0, xoffset=0.25, width=0.5, xrange=None):
        results = self.results
        a00 = results['PTC_A00']
        error = results['PTC_A00_ERROR']
        ymin = min(a00 - error*2)*1e6
        ymax = max(a00 + error*2)*1e6
        yname = 'Brighter-Fatter A00 (1e-6/e-)'
        win = plot.xyplot(results['AMP'], results['PTC_A00']*1e6,
                          yerr=results['PTC_A00_ERROR']*1e6, xname='AMP',
                          yname=yname, xrange=xrange, yrange=(ymin, ymax))
        win.set_title("PTC Brighter-Fatter, %s" % self.plot_title)

    def ptc_turnoff(self, oplot=0, xoffset=0.25, width=0.5, xrange=None):
        results = self.results
        yname = 'PTC Turnoff (DN)'
        win = plot.xyplot(results['AMP'], results['PTC_TURNOFF'],
                          xname='AMP', yname=yname, xrange=xrange)
        win.set_title("PTC Turnoff, %s" % self.plot_title)

    def noise(self, oplot=0, xoffset=0.2, width=0.2, color='b'):
        results = self.results
        read_noise = results['READ_NOISE']
        try:
            system_noise = results['SYSTEM_NOISE']
            total_noise = results['TOTAL_NOISE']
        except KeyError:
            system_noise = np.zeros(len(read_noise))
            total_noise = read_noise
        ymax = max(1.2*max(np.concatenate((read_noise,
                                           system_noise,
                                           total_noise))), 10)
        win = plot.bar(results['AMP'] - xoffset - xoffset/2., read_noise,
                       xname='Amp',
                       yname='Noise (rms e-)',
                       xrange=(0, len(read_noise)+1),
                       yrange=(0, ymax),
                       color=color, width=width)
        plot.bar(results['AMP'] - xoffset/2., system_noise,
                 oplot=1, color='g', width=width)
        plot.bar(results['AMP'] + xoffset - xoffset/2., total_noise,
                 oplot=1, color='c', width=width)
        plot.legend('bgc', ('Read Noise', 'System Noise', 'Total Noise'))
        plot.hline(8)
        win.set_title("Read Noise, %s" % self.plot_title)

    def total_noise(self, exptime=16, spec=9, dark95s=None):
        """
        Plot the total noise (electronic + 95th percentile dark current
        shot noise).
        """
        amp = self.results['AMP']
        total_noise = self.results['TOTAL_NOISE']
        if dark95s is None:
            shot_noise = self.results['DARK_CURRENT_95']*exptime
        else:
            shot_noise = np.array([dark95s[x] for x in amp])*exptime
        electronic_noise = np.sqrt(total_noise**2 + shot_noise)
        color_cycler = plt.rcParams['axes.prop_cycle']()
        npts = len(total_noise)
        dx = 0.075
        for noise, label, xoffset \
                in zip((total_noise, shot_noise, electronic_noise),
                       ('READ_NOISE', 'DC95_SHOT_NOISE', 'TOTAL_NOISE'),
                       (-dx*np.ones(npts), np.zeros(npts), dx*np.ones(npts))):
            color = next(color_cycler)['color']
            plt.plot(amp + xoffset, noise, '.', color=color, label=label)
        plt.xlabel('Amp')
        plt.ylabel('Noise (rms e-)')
        plt.title('Noise, %s' % self.plot_title)
        plt.legend()
        plt.plot([0, 17], [spec, spec], 'k:')
        axis = list(plt.axis())
        axis[1] = 17
        axis[-1] = max(10, axis[-1])
        plt.axis(axis)

    def full_well(self, gain_range=(1, 6), figsize=(11, 8.5),
                  ptc_file=None, detresp_file=None):
        if detresp_file is None:
            detresp_file = self._fullpath('%s_det_response.fits'
                                          % self.sensor_id)
        detresp = DetectorResponse(detresp_file, gain_range=gain_range)
        for amp in imutils.allAmps(detresp_file):
            subplot = self.subplot(amp)
            if amp == 1:
                win = plot.Window(subplot=subplot, figsize=figsize,
                                  xlabel=r'pd current $\times$ exposure',
                                  ylabel=r'$10^3$ e- per pixel', size='large')
                win.frameAxes.text(0.5, 1.08, 'Full Well, %s' % self.plot_title,
                                   horizontalalignment='center',
                                   verticalalignment='top',
                                   transform=win.frameAxes.transAxes,
                                   size='large')
            else:
                win.select_subplot(*subplot)
            self._offset_subplot(win)
            try:
                detresp.full_well(amp, make_plot=True, plotter=plot,
                                  multipanel=True)
            except (ValueError, RuntimeError, TypeError):
                pass
            pylab.annotate('Amp %i' % amp, (0.1, 0.8),
                           xycoords='axes fraction', size='x-small')

    @property
    def linearity_results(self, gain_range=(1, 6),
                          ptc_file=None, detresp_file=None):
        gain_range = self._gain_range
        ptc_file = self._ptc_file
        detresp_file = self._detresp_file
        if self._linearity_results is not None and detresp_file is None:
            return self._linearity_results
        self._linearity_results = {}
        if ptc_file is not None:
            ptc = fits.open(ptc_file)
        else:
            try:
                ptc = fits.open(self._fullpath('%s_ptc.fits' % self.sensor_id))
            except IOError:
                ptc = None
        if detresp_file is None:
            detresp_file = self._fullpath('%s_det_response.fits'
                                          % self.sensor_id)
        detresp = DetectorResponse(detresp_file, ptc=ptc, gain_range=gain_range)
        for amp in imutils.allAmps(detresp_file):
            try:
                self._linearity_results[amp] \
                    = detresp.linearity(amp, fit_range=self._Ne_bounds)
            except Exception as eObj:
                print("EOTestPlots.linearity: amp %i" % amp)
                print("  ", eObj)
        if ptc is not None:
            ptc.close()
        return self._linearity_results

    def linearity(self, gain_range=(1, 6), max_dev=0.02, figsize=(11, 8.5),
                  ptc_file=None, detresp_file=None, use_exptime=False,
                  Ne_bounds=(1e3, 9e4)):
        self._gain_range = gain_range
        self._ptc_file = ptc_file
        self._detresp_file = detresp_file
        self._Ne_bounds = Ne_bounds
        for amp in imutils.allAmps(detresp_file):
            #
            # Set up the plotting subwindow.
            subplot = self.subplot(amp)
            if use_exptime:
                xlabel = 'exposure time (s)'
            else:
                xlabel = r'pd current $\times$ exposure'
            if amp == 1:
                win = plot.Window(subplot=subplot, figsize=figsize,
                                  xlabel=xlabel, ylabel='', size='large')
                win.frameAxes.text(0.5, 1.08, 'Linearity, %s' % self.plot_title,
                                   horizontalalignment='center',
                                   verticalalignment='top',
                                   transform=win.frameAxes.transAxes,
                                   size='large')
            else:
                win.select_subplot(*subplot)

            #
            # Get the linearity fit for this amp.
            try:
                maxdev, fit_pars, Ne, flux, Ne_f, flux_f \
                    = self.linearity_results[amp][:6]
            except KeyError:
                continue
            #
            # Compute the fractional residuals
            f1 = np.poly1d(fit_pars)
            dNfrac = 1 - Ne/f1(flux)
            self._offset_subplot(win)
            # Resize subplot for plotting e-/pixel vs flux
            bbox = win.axes[-1].get_position()
            top_pts = bbox.get_points()
            bot_pts = copy.deepcopy(top_pts)
            dx, dy = top_pts[1] - top_pts[0]
            top_pts[0][1] += dy/4
            bbox.set_points(top_pts)
            win.axes[-1].set_position(bbox)
            #
            # Plot Ne vs flux
            try:
                win.axes[-1].loglog(flux, Ne, 'ko', markersize=3)
            except Exception as eObj:
                print("EOTestPlots.linearity: amp %i" % amp)
                print("  ", eObj)
            # Plot the fitted points
            try:
                win.axes[-1].loglog(flux_f, Ne_f, 'go', markersize=3)
            except Exception as eObj:
                print("EOTestPlots.linearity: amp %i" % amp)
                print("  ", eObj)
            try:
                win.axes[-1].loglog(flux, f1(flux), 'r-')
            except Exception as eObj:
                print("EOTestPlots.linearity: amp %i" % amp)
                print("  ", eObj)
            sys.stdout.flush()

            # Plot horizontal lines showing the range of the linearity
            # spec in e-/pixel.
            xmin, xmax, ymin, ymax = pylab.axis()
            win.axes[-1].loglog([xmin, xmax], [Ne_bounds[0], Ne_bounds[0]],
                                'k:')
            win.axes[-1].loglog([xmin, xmax], [Ne_bounds[1], Ne_bounds[1]],
                                'k:')

            # Label plots by amplifier number.
            pylab.annotate('Amp %i' % amp, (0.2, 0.8),
                           xycoords='axes fraction', size='x-small')
            if amp in (1, 5, 9, 13):
                win.axes[-1].set_ylabel('e-/pixel')
            for label in win.axes[-1].get_xticklabels():
                label.set_visible(False)

            # Add fractional residuals sub-subplot.
            bot_rect = [bot_pts[0][0], bot_pts[0][1], dx, dy/4.]
            bot_ax = win.fig.add_axes(bot_rect, sharex=win.axes[-1])
            bot_ax.semilogx(flux, dNfrac, 'ko', markersize=3)
            bot_ax.semilogx(flux, dNfrac, 'k:')
            bot_ax.semilogx(flux, np.zeros(len(Ne)), 'r-')
            pylab.locator_params(axis='y', nbins=5, tight=True)
            plot.setAxis(yrange=(-1.5*max_dev, 1.5*max_dev))

    def linearity_resids(self, gain_range=(1, 6), max_dev=0.02,
                         figsize=(11, 8.5), ptc_file=None, detresp_file=None,
                         Ne_bounds=(1e3, 9e4), use_exptime=False):
        self._gain_range = gain_range
        self._ptc_file = ptc_file
        self._detresp_file = detresp_file
        self._Ne_bounds = Ne_bounds
        for amp in imutils.allAmps(detresp_file):
            subplot = self.subplot(amp)
            if use_exptime:
                xlabel = 'exposure time (s)'
            else:
                xlabel = r'pd current $\times$ exposure'
            if amp == 1:
                win = plot.Window(subplot=subplot, figsize=figsize,
                                  xlabel=xlabel, ylabel='', size='large')
                win.frameAxes.text(0.5, 1.08,
                                   'Linearity residuals, %s' % self.plot_title,
                                   horizontalalignment='center',
                                   verticalalignment='top',
                                   transform=win.frameAxes.transAxes,
                                   size='large')
            else:
                win.select_subplot(*subplot)

            try:
                maxdev, fit_pars, Ne, flux = self.linearity_results[amp][:4]
            except KeyError:
                continue

            f1 = np.poly1d(fit_pars)
            dNfrac = 1 - Ne/f1(flux)

            self._offset_subplot(win)
            win.axes[-1].semilogx(flux, dNfrac, 'ko', markersize=3)
            win.axes[-1].semilogx(flux, dNfrac, 'k:')
            plot.setAxis(yrange=(-1.5*max_dev, 1.5*max_dev))
            xmin, xmax, ymin, ymax = pylab.axis()
            win.axes[-1].semilogx([xmin, xmax], [0, 0], 'r-')
            pylab.locator_params(axis='y', nbins=5, tight=True)
            # Plot max_dev range as horizontal lines.
            win.axes[-1].semilogx([xmin, xmax], [max_dev, max_dev], 'k--')
            win.axes[-1].semilogx([xmin, xmax], [-max_dev, -max_dev], 'k--')
            pylab.yticks([-max_dev, 0, max_dev])
            # Plot flux bounds corresponding to range of Ne values over which
            # the linearity spec is written
            flux_min = (Ne_bounds[0] - fit_pars[1])/fit_pars[0]
            flux_max = (Ne_bounds[1] - fit_pars[1])/fit_pars[0]
            print('amp, flux bounds, fit_pars:', amp, flux_min, flux_max, \
                fit_pars)
            win.axes[-1].semilogx([flux_min, flux_min], [ymin, ymax], 'k--')
            win.axes[-1].semilogx([flux_max, flux_max], [ymin, ymax], 'k--')
            # Label plots by amplifier number.
            pylab.annotate('Amp %i' % amp, (0.2, 0.9),
                           xycoords='axes fraction', size='x-small')
            if amp in (1, 5, 9, 13):
                # Add y-axis label for plots at left column.
                win.axes[-1].set_ylabel('e-/pixel residuals')
            else:
                # Suppress tick labels for other y-axes.
                for label in win.axes[-1].get_yticklabels():
                    label.set_visible(False)

    def cte_profiles(self, flux_level, sflat_file, mask_files,
                     figsize=(11, 8.5), serial=True):
        """
        Plot an array of serial or parallel cte profiles.
        """
        ccd = MaskedCCD(sflat_file, mask_files=mask_files)
        gains = dict((amp, gain) for amp, gain
                     in zip(self.results['AMP'], self.results['GAIN']))
        cti = dict((amp, Estimator()) for amp in ccd)
        bias_est = {}
        if serial:
            direction = 'Serial'
            xlabel = 'column #'
        else:
            direction = 'Parallel'
            xlabel = 'row #'

        keyname = '_'.join(('CTI', flux_level.upper(), direction.upper()))
        for amp in ccd:
            cti[amp].value = self.results[keyname][amp-1]
            cti[amp].error = self.results[keyname + '_ERROR'][amp-1]
            bias_est[amp] = gains[amp]*bias_estimate(ccd[amp], ccd.amp_geom,
                                                     serial=serial)
        title = '%s CTE profiles, %s flux, %s' % (direction, flux_level,
                                                  self.plot_title)
        for amp in ccd:
            subplot = self.subplot(amp)
            if amp == 1:
                win = plot.Window(subplot=subplot, figsize=figsize,
                                  xlabel=xlabel,
                                  ylabel='mean signal - bias (e-/pixel)',
                                  size='large')
                win.frameAxes.text(0.5, 1.08, title,
                                   horizontalalignment='center',
                                   verticalalignment='top',
                                   transform=win.frameAxes.transAxes,
                                   size='large')
            else:
                win.select_subplot(*subplot)
            self._offset_subplot(win)
            cte_profile(win.axes[-1], ccd[amp], gains[amp],
                        ccd.amp_geom, cti[amp], bias_est[amp], serial=serial)
            pylab.annotate('Amp %i\nCTI=%.2e\n    +/-%.2e'
                           % (amp, cti[amp].value, cti[amp].error),
                           (0.5, 0.75), xycoords='axes fraction',
                           size='x-small')

    def qe_ratio(self, ref, amp=None, qe_file=None):
        if qe_file is not None:
            self._qe_file = qe_file
        if amp is None:
            amps = imutils.allAmps(self._qe_file)
        else:
            amps = (amp,)
        for amp in amps:
            print("Amp", amp)
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
                             list(self.band_pass.keys()) if b in bands])
        band_wls_errs = np.array([(self.band_pass[b][1]-self.band_pass[b][0])/2.
                                  for b in list(self.band_pass.keys())
                                  if b in bands])
        wl = qe_data[1].data.field('WAVELENGTH')
        qe = {}
        for amp in imutils.allAmps(self._qe_file):
            qe[amp] = qe_data[1].data.field('AMP%02i' % amp)
            win = plot.curve(wl, qe[amp], xname='wavelength (nm)',
                             yname='QE (% e-/photon)', oplot=amp-1,
                             xrange=(300, 1100), yrange=(0, 120))
            if amp == 1:
                win.set_title('QE, %s' % self.plot_title)
            qe_band = qe_data[2].data.field('AMP%02i' % amp)
            plot.xyplot(band_wls, qe_band, xerr=band_wls_errs,
                        oplot=1, color='g')
        plot.hline(100)

    def flat_fields(self, lambda_dir, nsig=3, cmap=pylab.cm.hot, annotation=''):
        glob_string = os.path.join(lambda_dir, '*_lambda_*.fits')
        #print glob_string
        flats = sorted(glob.glob(glob_string))
        flats = [x for x in flats if x.find('bias') == -1]
        wls = []
        for flat in flats:
            try:
                wl = int(float(os.path.basename(flat).split('_')[2]))
            except ValueError:
                wl = int(float(os.path.basename(flat).split('_')[3]))
            wls.append(wl)
        wls = np.array(wls)
        #print wls
        # Loop over PRNU wavelengths and generate a png for each.
        gains = dict([(amp, gain) for amp, gain
                      in zip(self.results['AMP'], self.results['GAIN'])])
        for wl in self.prnu_wls:
            try:
                target = np.where(wls == wl)[0][0]
                win = plot_flat(flats[target], nsig=nsig, cmap=cmap, wl=wl,
                                gains=gains)
                plt.annotate(annotation, (1, -0.1), xycoords='axes fraction',
                             horizontalalignment='right',
                             verticalalignment='bottom')
                pylab.savefig('%s_%04inm_flat.png' % (self.sensor_id, wl))
                pylab.close()
            except IndexError:
                pass

    def confluence_tables(self, outfile=False):
        if outfile:
            output = open(self._outputpath('%s_results.txt'
                                           % self.sensor_id), 'w')
        else:
            output = sys.stdout
        # Write the per amp results.
        for name in self.results.colnames:
            output.write('|| %s' % name)
        output.write('||\n')
#        format = '| %i | %.2f | %.2f | %i | %.1e | %.1e | %.1e | %i | %i | %.1e | %.2f |\n'
        format = '| %i |' + (len(self.results.colnames)-1)*' %.2e |' + '\n'
        for i, amp in enumerate(self.results['AMP']):
            output.write(format % tuple([self.results[x][i]
                                         for x in self.results.colnames]))
        output.write('\n')
        # Write the CCD-wide results.
        # PRNU:
        prnu_results = self.prnu_results
        output.write("|| wavelength || stdev of pixel values || mean || stdev/mean ||\n")
        for wl, stdev, mean in zip(prnu_results['WAVELENGTH'],
                                   prnu_results['STDEV'], prnu_results['MEAN']):
            if stdev > 0:
                output.write("| %i | %12.4e | %12.4e | %12.4e |\n"
                             % (wl, stdev, mean, stdev/mean))
            else:
                output.write("| %i | ... | ... | ... |\n" % wl)
        if outfile:
            output.close()

    @property
    def prnu_results(self):
        my_prnu_results = fits.open(self.results.infile)['PRNU_RESULTS'].data
        return my_prnu_results

    def latex_table(self, outfile=None, hspace=None):
        lines = []
        lines.append(self.specs.latex_header(hspace=hspace))
        for spec in self.specs:
            lines.append(self.specs[spec].latex_entry())
        lines.append(self.specs.latex_footer())
        my_table = '\n'.join(lines) + '\n'
        if outfile is not None:
            output = open(outfile, 'w')
            output.write(my_table)
            output.close()
        return my_table


class CcdSpecs(OrderedDict):
    _job_name_map = dict(read_noise=('CCD-007',),
                         read_noise_offline=('CCD-007',),
                         flat_pairs=('CCD-008', 'CCD-009'),
                         flat_pairs_offline=('CCD-008', 'CCD-009'),
                         cte=('CCD-010', 'CCD-011'),
                         cte_offline=('CCD-010', 'CCD-011'),
                         bright_defects=('CCD-012a', 'CCD-012c'),
                         bright_defects_offline=('CCD-012a', 'CCD-012c'),
                         dark_defects=('CCD-012b', 'CCD-012d'),
                         dark_defects_offline=('CCD-012b', 'CCD-012d'),
                         traps=('CCD-012e',),
                         traps_offline=('CCD-012e',),
                         dark_current=('CCD-014',),
                         dark_current_offline=('CCD-014',),
                         qe_analysis=('CCD-021', 'CCD-022', 'CCD-023',
                                      'CCD-024', 'CCD-025', 'CCD-026'),
                         qe_offline=('CCD-021', 'CCD-022', 'CCD-023',
                                     'CCD-024', 'CCD-025', 'CCD-026'),
                         prnu=('CCD-027',),
                         prnu_offline=('CCD-027',),
                         fe55_analysis=('CCD-028',),
                         fe55_offline=('CCD-028',),
                         )

    def __init__(self, results_file, xtalk_file=None, plotter=None,
                 prnu_wls=()):
        super(CcdSpecs, self).__init__()
        self.plotter = plotter
        self.prnu_wls = prnu_wls
        self.prnu_specs = OrderedDict()
        self._createSpecs()
        try:
            self._ingestResults(results_file, xtalk_file=xtalk_file)
        except Exception as eobj:
            pass
#            print("EOTestPlots.CcdSpecs: exception:")
#            print("  ", str(eobj))

    def add_job_ids(self, summary_files):
        for summary_lims_file in summary_files:
            foo = json.loads(open(summary_lims_file).read())
            for result in foo:
                if 'job_id' in result:
                    try:
                        specids = self._job_name_map[result['job_name']]
                        for specid in specids:
                            self[specid].job_id = result['job_id']
                            if specid == 'CCD-027':
                                for prnu_spec in list(self.prnu_specs.values()):
                                    prnu_spec.job_id = result['job_id']
                    except KeyError:
                        pass

    def factory(self, *args, **kwds):
        spec = CcdSpec(*args, **kwds)
        self[spec.name] = spec
        return spec

    def _createSpecs(self):
        self.factory('CCD-007', 'Read Noise', spec='$< 8$\,\electron rms')
        self.factory('CCD-008', 'Blooming Full Well',
                     spec='$<175000$\,\electron')
        self.factory('CCD-009', 'Nonlinearity', spec='$<2\\%$')
        self.factory('CCD-010', 'Serial CTE', spec='$> 1 - \\num{5e-6}$')
        self.factory('CCD-011', 'Parallel CTE', spec='$> 1 - \\num{3e-6}$')
        self.factory('CCD-012', '\\twolinecell{Active Imaging Area \\\\and Cosmetic Quality}',
                     spec='\\twolinecell{$<0.5$\\% defective \\\\pixels}')
        self.factory('CCD-012a', 'Bright Pixels')
        self.factory('CCD-012b', 'Dark Pixels')
        self.factory('CCD-012c', 'Bright Columns')
        self.factory('CCD-012d', 'Dark Columns')
        self.factory('CCD-012e', 'Traps')
        self.factory('CCD-013', 'Crosstalk', spec='$<0.19$\\%')
        self.factory('CCD-014',
                     '\\twolinecell{Dark Current \\\\95th Percentile}',
                     spec='$<0.2$\,\electron\,s$^{-1}$')
        self.factory('CCD-021', 'u Band QE', spec='$> 41$\\%')
        self.factory('CCD-022', 'g Band QE', spec='$> 78$\\%')
        self.factory('CCD-023', 'r Band QE', spec='$> 83$\\%')
        self.factory('CCD-024', 'i Band QE', spec='$> 82$\\%')
        self.factory('CCD-025', 'z Band QE', spec='$> 75$\\%')
        self.factory('CCD-026', 'y Band QE', spec='$> 21$\\%')
        self.factory('CCD-027', 'PRNU', spec='$<5$\\%')
        for wl in self.prnu_wls:
            self.prnu_specs[wl] = CcdSpec("CCD-027 (%inm)" % wl, 'PRNU',
                                          spec='$<5$\\%')
        self.factory('CCD-028', 'Point Spread Function', spec='$\sigma < 5\mu$')

    @staticmethod
    def latex_header(hspace=None):
        return CcdSpec.latex_header(hspace=hspace)

    @staticmethod
    def latex_footer():
        return CcdSpec.latex_footer()

    def latex_table(self, hspace=None):
        output = []
        output.append(self.latex_header(hspace=hspace))
        for name, spec in list(self.items()):
            output.append(spec.latex_entry())
        output.append(self.latex_footer())
        return '\n'.join(output) + '\n'

    def _ingestResults(self, results_file, xtalk_file=None):
        self.results = EOTestResults(results_file)
        rn = self.results['READ_NOISE']
        self['CCD-007'].measurement = '$%s$--$%s$\,\electron\,rms' \
            % (min_str(rn, '%.2f'), max_str(rn, '%.2f'))
        self['CCD-007'].ok = (max(rn) < 8)
        fw = self.results['FULL_WELL']
        self['CCD-008'].measurement = '$%s$--$%s$\,\electron' % (min_str(fw, '%i'), max_str(fw, '%i'))
        self['CCD-008'].ok = (max(fw) < 175000)
        max_frac_dev = self.results['MAX_FRAC_DEV']
        self['CCD-009'].measurement = '\\twolinecell{max. fractional deviation \\\\from linearity: $\\num{%s}$}' % max_str(
            max_frac_dev, '%.1e')
        self['CCD-009'].ok = (max(max_frac_dev) < 0.02)
        scti = self.results['CTI_HIGH_SERIAL']
        scti = np.concatenate((scti, self.results['CTI_LOW_SERIAL']))
        try:
            self.results['CTI_HIGH_SERIAL_ERROR']
            scti_err = self.results['CTI_HIGH_SERIAL_ERROR']
            scti_err = np.concatenate((scti_err, self.results['CTI_LOW_SERIAL_ERROR']))
        except KeyError:
            scti_err = np.zeros(32)
        self['CCD-010'].measurement = '$1%s$ (min. value)' % latex_minus_max(scti, scti_err)
        self['CCD-010'].ok = (max(scti) < 5e-6)
        pcti = self.results['CTI_HIGH_PARALLEL']
        pcti = np.concatenate((pcti, self.results['CTI_LOW_PARALLEL']))
        try:
            pcti_err = self.results['CTI_HIGH_PARALLEL_ERROR']
            pcti_err = np.concatenate((pcti_err, self.results['CTI_LOW_PARALLEL_ERROR']))
        except KeyError:
            pcti_err = np.zeros(32)
        self['CCD-011'].measurement = '$1%s$ (min. value)' % latex_minus_max(pcti, pcti_err)
        self['CCD-011'].ok = (max(pcti) < 3e-6)
        num_bright_pixels = sum(self.results['NUM_BRIGHT_PIXELS'])
        self['CCD-012a'].measurement = '%i' % num_bright_pixels
        num_dark_pixels = sum(self.results['NUM_DARK_PIXELS'])
        self['CCD-012b'].measurement = '%i' % num_dark_pixels
        num_bright_columns = sum(self.results['NUM_BRIGHT_COLUMNS'])
        self['CCD-012c'].measurement = '%i' % num_bright_columns
        num_dark_columns = sum(self.results['NUM_DARK_COLUMNS'])
        self['CCD-012d'].measurement = '%i' % num_dark_columns
        num_traps = sum(self.results['NUM_TRAPS'])
        self['CCD-012e'].measurement = '%i' % num_traps

        try:
            num_pixels = (self.results['TOTAL_NUM_PIXELS']
                          - self.results['ROLLOFF_MASK_PIXELS'])
        except Exception:
            num_pixels = 16129000
        col_size = 2002 - 9 # exclude masked edge rolloff.
        num_defects = (num_bright_pixels + num_dark_pixels + num_traps
                       + col_size*(num_bright_columns + num_dark_columns))
        self['CCD-012'].measurement = 'defective pixels: %i (%.4f\\%%)' \
            % (num_defects, 100.*float(num_defects)/float(num_pixels))
        self['CCD-012'].ok = (float(num_defects)/float(num_pixels) < 5e-3)
        if xtalk_file is not None:
            ctm = CrosstalkMatrix(xtalk_file)
            for i in range(len(ctm.matrix)):
                ctm.matrix[i][i] = 0
            crosstalk = max(np.abs(ctm.matrix.flatten()))
            self['CCD-013'].measurement = 'max. value: $\\num{%.2e}$\\%%' % (crosstalk*100.)
            self['CCD-013'].ok = (crosstalk < 1.9e-3)
        try:
            dark_current = self.results.output['AMPLIFIER_RESULTS'].header['DARK95']
        except KeyError:
            dark_current = max(self.results['DARK_CURRENT_95'])
        self['CCD-014'].measurement = '$\\num{%.2e}$\electron\,s$^{-1}$' % dark_current
        self['CCD-014'].ok = (dark_current < 0.2)

        bands = self.plotter.qe_data['QE_BANDS'].data.field('BAND')
        bands = OrderedDict([(band, []) for band in bands])
        for amp in self.results['AMP']:
            try:
                values = self.plotter.qe_data['QE_BANDS'].data.field('AMP%02i' % amp)
                for band, value in zip(bands, values):
                    bands[band].append(value)
            except KeyError:
                pass
        for band, specnum, minQE in zip('ugrizy', list(range(21, 27)),
                                        (41, 78, 83, 82, 75, 21)):
            try:
                qe_mean = np.mean(bands[band])
                self['CCD-0%i' % specnum].measurement = '%.1f\\%%' % qe_mean
                self['CCD-0%i' % specnum].ok = (qe_mean > minQE)
            except KeyError:
                self['CCD-0%i' % specnum].measurement = 'No data'
        prnu_results = self.plotter.prnu_results
        target_wls = set(EOTestPlots.prnu_wls)
        ratios = {}
        for wl, stdev, mean in zip(prnu_results['WAVELENGTH'],
                                   prnu_results['STDEV'], prnu_results['MEAN']):
            if stdev > 0:
                target_wls.remove(int(wl))
                ratios[wl] = stdev/mean
                if wl in self.prnu_specs:
                    if ratios[wl] < 0.01:
                        self.prnu_specs[wl].measurement = \
                            "\\num{%.1e}\\%%" % (ratios[wl]*100)
                    else:
                        self.prnu_specs[wl].measurement = \
                            "%.2f\\%%" % (ratios[wl]*100)
                    self.prnu_specs[wl].ok = (ratios[wl] < 5e-2)
        max_ratio = max(ratios.values())
        max_wl = list(ratios.keys())[np.where(list(ratios.values()) == max_ratio)[0][0]]
        if max_ratio < 0.01:
            self['CCD-027'].measurement = 'max. variation = \\num{%.1e}\\%% at %i\\,nm' % (
                max_ratio*100, max_wl)
        else:
            self['CCD-027'].measurement = 'max. variation = %.2f\\%% at %i\\,nm' % (max_ratio*100, max_wl)
        if target_wls:
            measurement = self['CCD-027'].measurement + '\\\\missing wavelengths: ' \
                + ', '.join([str(x) for x in sorted(target_wls)]) + "\\,nm"
            self['CCD-027'].measurement = '\\twolinecell{%s}' % measurement
        self['CCD-027'].ok = (max_ratio < 5e-2)
        psf_sigma = max(self.results['PSF_SIGMA'])
        self['CCD-028'].measurement = '$%.2f\,\mu$ (max. value)' % psf_sigma
        self['CCD-028'].ok = (psf_sigma < 5.)


class CcdSpec(object):
    _latex_status = dict([(True, '\ok'), (False, '\\fail'), (None, '$\cdots$')])

    def __init__(self, name, description, spec=None, ok=None, measurement=None):
        self.name = name
        self.description = description
        self.spec = spec
        self.ok = ok
        self.measurement = measurement
        self.job_id = '$\cdots$'

    @staticmethod
    def _table_cell(value):
        if value is None:
            return '$\cdots$'
        else:
            return value

    @staticmethod
    def latex_header(hspace=None):
        if hspace is None:
            header = """\\begin{table}[h]
\\centering
\\begin{tabular}{|c|l|l|l|l|l|}
\hline
\\textbf{Status} & \\textbf{Spec. ID} & \\textbf{Description} & \\textbf{Specification} & \\textbf{Measurement}  & \\textbf{Job ID} \\\\ \hline"""
        else:
            header = """\\begin{table}[h]
\\hspace{%s}
\\begin{tabular}{|c|l|l|l|l|l|}
\hline
\\textbf{Status} & \\textbf{Spec. ID} & \\textbf{Description} & \\textbf{Specification} & \\textbf{Measurement}  & \\textbf{Job ID} \\\\ \hline""" % hspace
        return header

    @staticmethod
    def latex_footer():
        footer = "\end{tabular}\n\end{table}"
        return footer

    def latex_entry(self):
        entry = '%s & %s & %s & %s & %s & %s \\\\ \hline' % \
                (self._latex_status[self.ok],
                 self.name,
                 self.description,
                 self._table_cell(self.spec),
                 self._table_cell(self.measurement),
                 self.job_id)
        return entry

    def latex_table(self):
        lines = []
        lines.append(CcdSpecs.latex_header())
        lines.append(self.latex_entry())
        lines.append(CcdSpecs.latex_footer())
        return '\n'.join(lines)
