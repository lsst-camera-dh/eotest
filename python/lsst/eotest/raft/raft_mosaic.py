"""
Code to perform raft-level mosaicking from single sensor frames
compliant with LCA-13501.
"""
import os
import copy
from collections import defaultdict
import numpy as np
import astropy.io.fits as fits
import matplotlib
import matplotlib.pyplot as plt
import lsst.eotest.image_utils as imutils
import lsst.eotest.sensor as sensorTest
from lsst.eotest.sensor.EOTestPlots import cmap_range


__all__ = ['RaftMosaic', 'CornerRaftMosaic', 'make_raft_mosaic']


def make_raft_mosaic(fits_files, gains=None, bias_subtract=True,
                     segment_processor=None, bias_frames=None,
                     dark_currents=None):
    corner_raft_slots = {'SW0', 'SW1', 'SG0', 'SG1'}
    if corner_raft_slots.intersection(fits_files):
        return CornerRaftMosaic(fits_files, bias_frames=bias_frames)
    return RaftMosaic(fits_files, gains=gains,
                      bias_subtract=bias_subtract,
                      segment_processor=segment_processor,
                      bias_frames=bias_frames,
                      dark_currents=dark_currents)


class RaftMosaic:
    """
    Raft-level mosaic of individual CCDs in Science Rafts.
    """

    def __init__(self, fits_files, gains=None, bias_subtract=True,
                 nx=12700, ny=12700, nx_segments=8, ny_segments=2,
                 segment_processor=None, bias_frames=None,
                 dark_currents=None, e2v_xoffset=21):
        """
        Parameters
        ----------
        fits_files : dict
            Dictionary of single sensor FITS files, keyed by raft slot
            name.  These files should conform to LCA-13501.
        gains : dict [None]
            Dictionary (keyed by slot name) of dictionaries (one per
            FITS file) of system gain values for each amp.  If
            None, then do not apply gain correction.
        bias_subtract : bool [True]
            Flag do to a bias subtraction based on the serial overscan
            or provided bias frame.
        nx : int [12700]
            Number of pixels in the x (serial) direction.
        ny : int, [12700]
            Number of pixels in the y (parallel) direction.
        nx_segments : int [8]
            Number of segments in the x (serial) direction.
        ny_segments : int [2]
            Number of pixels in the y (parallel) direction.
        segment_processor : function [None]
            Function to apply to pixel data in each segment. If None,
            then set do the standard bias subtraction and gain correction.
        bias_frames : dict [None]
            Dictionary of single sensor bias frames, keyed by raft slot.
            If None, then just do the bias level subtraction using
            overscan region.
        dark_currents : dict [None]
            Dictionary of dictionaries of dark current values per amp
            in e-/s, keyed by raft slot and by amp number. If None, then
            dark current subtraction is not applied.
        e2v_xoffset : int [21]
            Offset in serial direction for CRVAL1Q parameter to get
            07-00 segments to align properly with 10-17 segments
            in BOT-level FITS files for e2V CCDs.
        """
        self.fits_files = fits_files
        with fits.open(list(fits_files.values())[0]) as hdu_list:
            self.raft_name = hdu_list[0].header['RAFTNAME']
            try:
                self.wl = hdu_list[0].header['MONOWL']
            except KeyError:
                self.wl = 0
        self.image_array = np.zeros((nx, ny), dtype=np.float32)
        self.nx = nx
        self.ny = ny
        self.nx_segments = nx_segments
        self.ny_segments = ny_segments
        self.segment_processor = segment_processor
        self._amp_coords = defaultdict(dict)
        if gains is None:
            # Assume unit gain for all amplifiers.
            unit_gains = dict([(i, 1) for i in range(1, 17)])
            gains = dict([(slot, unit_gains) for slot in fits_files])
        for slot, filename in list(fits_files.items()):
            #print("processing", os.path.basename(filename))
            bias_frame = bias_frames[slot] if bias_frames is not None else None
            ccd = sensorTest.MaskedCCD(filename, bias_frame=bias_frame)
            if dark_currents is not None:
                try:
                    dark_time = ccd.md.get('DARKTIME')
                except:
                    dark_time = ccd.md.get('EXPTIME')
            with fits.open(filename) as hdu_list:
                for amp, hdu in zip(ccd, hdu_list[1:]):
                    if dark_currents:
                        dark_correction = dark_time*dark_currents[slot][amp]
                    else:
                        dark_correction = 0
                    self._set_segment(slot, ccd, amp, hdu, gains[slot][amp],
                                      bias_subtract, dark_correction,
                                      e2v_xoffset=e2v_xoffset)

    def _set_segment(self, slot, ccd, amp, hdu, amp_gain, bias_subtract,
                     dark_correction, e2v_xoffset=21):
        """
        Set the pixel values in the mosaic from the segment values.
        """
        # Get the trimmed masked image, with or without bias subtraction.
        if bias_subtract:
            mi = ccd.unbiased_and_trimmed_image(amp)
        else:
            mi = ccd[amp].Factory(ccd[amp], ccd.amp_geom.imaging)
        # Apply gain correction.
        seg_array = np.array(amp_gain*copy.deepcopy(mi.getImage().getArray()),
                             dtype=np.float32)
        # Apply dark current correction.
        seg_array -= dark_correction
        # Determine flip in serial direction based on 1, 1 element of
        # transformation matrix.
        if hdu.header['PC1_1Q'] < 0:
            seg_array = seg_array[:, ::-1]
            xmax = int(hdu.header['CRVAL1Q'])
            xmin = xmax - ccd.amp_geom.nx
        else:
            xmin = int(hdu.header['CRVAL1Q'])
            if ccd.amp_geom.vendor == 'E2V':
                xmin += e2v_xoffset
            xmax = xmin + ccd.amp_geom.nx
        # Determine flip in parallel direction based on 2, 2 element
        # of transformation matrix.
        if hdu.header['PC2_2Q'] < 0:
            seg_array = seg_array[::-1, :]
            ymax = int(hdu.header['CRVAL2Q'])
            ymin = ymax - ccd.amp_geom.ny
        else:
            ymin = int(hdu.header['CRVAL2Q'])
            ymax = ymin + ccd.amp_geom.ny
        # Save coordinates of segment for later use.
        self._amp_coords[slot][amp] = xmin, xmax, ymin, ymax

        # Write the segment pixel values into the full raft image mosaic.
        if self.segment_processor is None:
            self.image_array[ymin:ymax, xmin:xmax] = seg_array
        else:
            xy_bounds = (xmin, xmax, ymin, ymax)
            self.image_array[ymin:ymax, xmin:xmax] = \
                self.segment_processor(slot, ccd, amp, xy_bounds=xy_bounds)

    def plot(self, title=None, cmap=plt.cm.hot, nsig=5, figsize=(10, 10),
             binsize=10, flipx=True, textcolor='c', annotation='',
             rotate180=False, vrange=None):
        """
        Render the raft mosaic.

        Parameters
        ----------
        title : str, optional
            The plot title. If None (default), then build the title
            from the RAFTNAME and MONOWL primary header keyword values.
        cmap : matplotlib.colors.Colormap, optional
            The color map to use. Default: matplotlib.pyplot.cm.hot.
        nsig : float, optional
            The n-sigma value for the sigma clipping used to determine
            the pixel value range over which the color map is mapped.
        figsize : (float, float), optional
            The width x height size of the figure in inches. Default: (10, 10).
        binsize : int, optional
            Rebin the plotted image data by binsize*binsize,
            averging over the coarser bin.  Default: 10
        flipx : bool, optional
            Flip full raft mosaic in x so that parity of image matches
            LCA-13381. Default: True
        textcolor : str, optional
            Color of the text for the segment and sensor labeling.
            Default: 'c' (cyan)
        annotation : str, optional
            Description of the plot, e.g., pixel units (ADU or e-),
            gain-corrected, bias-subtracted.  Default: ''
        rotate180 : bool [False]
            Flag to rotate the mosaic by 180 degrees to match the
            orientation of the focalplane mosiacs created for the
            BOT-level plots.
        vrange : (float, float) [None]
            Range of pixel values to plot.  If None, then the cmap_range
            function will be used.
        """
        plt.rcParams['figure.figsize'] = figsize
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        output_array = imutils.rebin_array(self.image_array, binsize,
                                           use_mean=True)
        if flipx:
            output_array = output_array[:, ::-1]
        if rotate180:
            ny, nx = output_array.shape
            rotated_array = np.zeros((nx, ny), dtype=output_array.dtype)
            for j in range(ny):
                rotated_array[:, ny-1-j] = output_array[::-1, j]
            output_array = rotated_array
        image = ax.imshow(output_array, interpolation='nearest', cmap=cmap)
        if vrange is None:
            # Set range and normalization of color map based on
            # sigma-clip of pixel values.
            vmin, vmax = cmap_range(output_array, nsig=nsig)
        else:
            vmin, vmax = vrange
        norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
        image.set_norm(norm)
        if title is None:
            title = "%s, %i nm" % (self.raft_name, self.wl)
        ax.set_title(title)
        fig.colorbar(image)
        # Turn off ticks and tick labels for x- and y-axes.
        plt.tick_params(axis='both', which='both',
                        top='off', bottom='off', left='off', right='off',
                        labelbottom='off', labelleft='off')
        # Label segments by sensor bay and segment number.
        for slot in self.fits_files:
            seg_coords = list(self._amp_coords[slot].values())[-8]
            xmin, xmax, ymin, ymax = seg_coords
            xx = float(xmax + xmin)/2./float(self.nx)
            if flipx:
                xx = 1 - xx
            yy = 1. - (float(ymax - ymin)*0.05 + ymin)/float(self.ny)
            if rotate180:
                xx = 1 - xx - 7*np.abs(xmax - xmin)/float(self.nx)
                yy = 1 - yy + 1.9*np.abs(ymax - ymin)/float(self.ny)
            plt.annotate('%s' % slot,
                         (xx, yy), xycoords='axes fraction',
                         size='x-small', horizontalalignment='center',
                         verticalalignment='center', color=textcolor)
            for amp, seg_coords in list(self._amp_coords[slot].items()):
                xmin, xmax, ymin, ymax = seg_coords
                xx = float(xmax + xmin)/2./float(self.nx)
                if flipx:
                    xx = 1. - xx
                if amp <= 8:
                    yy = 1. - (float(ymax - ymin)*0.85 + ymin)/float(self.ny)
                else:
                    yy = 1. - (float(ymax - ymin)*0.15 + ymin)/float(self.ny)
                if rotate180:
                    xx = 1 - xx
                    yy = 1 - yy
                plt.annotate('%s' % imutils.channelIds[amp],
                             (xx, yy), xycoords='axes fraction',
                             size='x-small', horizontalalignment='center',
                             verticalalignment='center', color=textcolor)
        plt.annotate(annotation, (1, -0.1), xycoords='axes fraction',
                     horizontalalignment='right', verticalalignment='bottom')
        return fig


class CornerRaftMosaic:
    # The amp_llc data are the nominal locations in pixel coordinates
    # of the lower left corner of each imaging segment for the
    # corresponding amp in the guide and wavefront sensors.  These
    # values assume a coordinate system which is rotated so that the
    # wavefront sensors are in the lower left corner of the raft.
    # They are based on sheet 4 of LCA-13381.
    amp_llc = dict()
    amp_llc['SG0'] = {amp: (2050, 4663 + (amp - 1)*509) for amp in range(1, 9)}
    amp_llc['SG0'].update({amp: (49, 8226 - (amp - 9)*509)
                           for amp in range(9, 17)})
    amp_llc['SG1'] = {amp: (8226 - (amp - 1)*509, 2050) for amp in range(1, 9)}
    amp_llc['SG1'].update({amp: (4663 + (amp - 9)*509, 49)
                           for amp in range(9, 17)})
    amp_llc['SW0'] = {amp: (3777 - (amp-1)*509, 88) for amp in range(1, 9)}
    amp_llc['SW1'] = {amp: (214 + (amp-1)*509, 2413) for amp in range(1, 9)}

    def __init__(self, fits_files, nx=11630, ny=11630, bias_frames=None):
        self.fits_files = fits_files
        with fits.open(list(fits_files.values())[0]) as hdus:
            self.raft_name = hdus[0].header['RAFTNAME']
            try:
                self.wl = hdus[0].header['MONOWL']
            except KeyError:
                self.wl = None
        self.image_array = np.zeros((nx, ny), dtype=np.float32)
        for slot, filename in fits_files.items():
            if slot not in self.amp_llc:
                continue
            bias_frame = bias_frames[slot] if bias_frames is not None else None
            ccd = sensorTest.MaskedCCD(filename, bias_frame=bias_frame)
            for amp in ccd:
                if slot.startswith('SW'):
                    self._set_sw_segment(slot, ccd, amp)
                else:
                    self._set_sg_segment(slot, ccd, amp)

    def _set_sw_segment(self, slot, ccd, amp):
        xmin, ymin = self.amp_llc[slot][amp]
        xmax = xmin + ccd.amp_geom.nx
        ymax = ymin + ccd.amp_geom.ny
        seg_array = np.array(copy.deepcopy(ccd.unbiased_and_trimmed_image(amp)
                                           .getImage().getArray()),
                             dtype=np.float32)
        if slot == 'SW0':
            # Flip amps in serial direction and parallel directions.
            # This is equivalent to a 180 degree rotation about the
            # segment center.
            seg_array = seg_array[::-1, ::-1]
        self.image_array[ymin:ymax, xmin:xmax] = seg_array

    def _set_sg_segment(self, slot, ccd, amp):
        xmin, ymin = self.amp_llc[slot][amp]
        if slot == 'SG0':
            # SGO sensors have their x- and y-coordinates swapped wrt
            # SG1 and wavefront and science raft sensors, so swap nx
            # and ny from the amp geometries.  The corresponding
            # transpose of the pixel data is performed below.
            xmax = xmin + ccd.amp_geom.ny
            ymax = ymin + ccd.amp_geom.nx
        else:
            xmax = xmin + ccd.amp_geom.nx
            ymax = ymin + ccd.amp_geom.ny

        seg_array = np.array(copy.deepcopy(ccd.unbiased_and_trimmed_image(amp)
                                           .getImage().getArray()),
                             dtype=np.float32)
        if amp < 9:
            seg_array = seg_array[::-1, :]
        if slot == 'SG0':
            seg_array = seg_array.transpose()
        else:
            seg_array = seg_array[:, ::-1]
        self.image_array[ymin:ymax, xmin:xmax] = seg_array

    def plot(self, title=None, cmap=plt.cm.hot, nsig=5, figsize=(10, 10),
             binsize=10, flipx=True, textcolor='c', annotation='',
             rotate180=False, vrange=None):
        """
        Render the raft mosaic.

        Parameters
        ----------
        title : str, optional
            The plot title. If None (default), then build the title
            from the RAFTNAME and MONOWL primary header keyword values.
        cmap : matplotlib.colors.Colormap, optional
            The color map to use. Default: matplotlib.pyplot.cm.hot.
        nsig : float, optional
            The n-sigma value for the sigma clipping used to determine
            the pixel value range over which the color map is mapped.
        figsize : (float, float), optional
            The width x height size of the figure in inches. Default: (10, 10).
        binsize : int, optional
            Rebin the plotted image data by binsize*binsize,
            averging over the coarser bin.  Default: 10
        flipx : bool, optional
            Flip full raft mosaic in x so that parity of image matches
            LCA-13381. Default: True
        textcolor : str, optional
            Color of the text for the segment and sensor labeling.
            Default: 'c' (cyan)
        annotation : str, optional
            Description of the plot, e.g., pixel units (ADU or e-),
            gain-corrected, bias-subtracted.  Default: ''
        rotate180 : bool [False]
            Flag to rotate the mosaic by 180 degrees to match the
            orientation of the focalplane mosiacs created for the
            BOT-level plots.
        vrange : (float, float) [None]
            Range of pixel values to plot.  If None, then the cmap_range
            function will be used.
        """
        plt.rcParams['figure.figsize'] = figsize
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        output_array = imutils.rebin_array(self.image_array, binsize,
                                           use_mean=True)
        if flipx:
            output_array = output_array[:, ::-1]
        if rotate180:
            ny, nx = output_array.shape
            rotated_array = np.zeros((nx, ny), dtype=output_array.dtype)
            for j in range(ny):
                rotated_array[:, ny-1-j] = output_array[::-1, j]
            output_array = rotated_array
        image = ax.imshow(output_array, interpolation='nearest', cmap=cmap)
        if vrange is None:
            # Set range and normalization of color map based on
            # sigma-clip of pixel values.
            vmin, vmax = cmap_range(output_array, nsig=nsig)
        else:
            vmin, vmax = vrange
        norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
        image.set_norm(norm)
        if title is None:
            if self.wl is None:
                title = self.raft_name
            else:
                title = "%s, %i nm" % (self.raft_name, self.wl)
        ax.set_title(title)
        fig.colorbar(image)
        # Turn off ticks and tick labels for x- and y-axes.
        plt.tick_params(axis='both', which='both',
                        top='off', bottom='off', left='off', right='off',
                        labelbottom='off', labelleft='off')
#        # Label segments by sensor bay and segment number.
#        for slot in self.fits_files:
#            seg_coords = list(self._amp_coords[slot].values())[-8]
#            xmin, xmax, ymin, ymax = seg_coords
#            xx = float(xmax + xmin)/2./float(self.nx)
#            if flipx:
#                xx = 1 - xx
#            yy = 1. - (float(ymax - ymin)*0.05 + ymin)/float(self.ny)
#            if rotate180:
#                xx = 1 - xx - 7*np.abs(xmax - xmin)/float(self.nx)
#                yy = 1 - yy + 1.9*np.abs(ymax - ymin)/float(self.ny)
#            plt.annotate('%s' % slot,
#                         (xx, yy), xycoords='axes fraction',
#                         size='x-small', horizontalalignment='center',
#                         verticalalignment='center', color=textcolor)
#            for amp, seg_coords in list(self._amp_coords[slot].items()):
#                xmin, xmax, ymin, ymax = seg_coords
#                xx = float(xmax + xmin)/2./float(self.nx)
#                if flipx:
#                    xx = 1. - xx
#                if amp <= 8:
#                    yy = 1. - (float(ymax - ymin)*0.85 + ymin)/float(self.ny)
#                else:
#                    yy = 1. - (float(ymax - ymin)*0.15 + ymin)/float(self.ny)
#                if rotate180:
#                    xx = 1 - xx
#                    yy = 1 - yy
#                plt.annotate('%s' % imutils.channelIds[amp],
#                             (xx, yy), xycoords='axes fraction',
#                             size='x-small', horizontalalignment='center',
#                             verticalalignment='center', color=textcolor)
#        plt.annotate(annotation, (1, -0.1), xycoords='axes fraction',
#                     horizontalalignment='right', verticalalignment='bottom')
        return fig
