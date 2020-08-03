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
                     dark_currents=None, nx=None, ny=None):
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
    nx : int [None]
        Number of pixels in the x (nominally serial) direction.  If
        None, then use defaults of RaftMosaic (12700) and
        CornerRaftMosaic (11630).
    ny : int [None]
        Number of pixels in the y (nominally parallel) direction.  If
        None, then use defaults of RaftMosaic (12700) and
        CornerRaftMosaic (11630).

    Returns
    -------
    RaftMosaic or CornerRaftMosaic
    """
    corner_raft_slots = {'SW0', 'SW1', 'SG0', 'SG1'}
    if corner_raft_slots.intersection(fits_files):
        return CornerRaftMosaic(fits_files, gains=gains,
                                bias_subtract=bias_subtract,
                                nx=nx, ny=ny,
                                bias_frames=bias_frames,
                                dark_currents=dark_currents)
    return RaftMosaic(fits_files, gains=gains,
                      bias_subtract=bias_subtract, nx=nx, ny=ny,
                      segment_processor=segment_processor,
                      bias_frames=bias_frames,
                      dark_currents=dark_currents)


class RaftMosaic:
    """
    Raft-level mosaic of individual CCDs in Science Rafts.
    """

    def __init__(self, fits_files, gains=None, bias_subtract=True,
                 nx=None, ny=None, segment_processor=None,
                 bias_frames=None, dark_currents=None, e2v_xoffset=21):
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
        nx : int [None]
            Number of pixels in the x (serial) direction. If None, use 12700.
        ny : int [None]
            Number of pixels in the y (parallel) direction. If None, use 12700.
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
                self.wl = None
        if nx is None:
            nx = 12700
        if ny is None:
            ny = 12700
        self.image_array = np.zeros((nx, ny), dtype=np.float32)
        self.nx = nx
        self.ny = ny
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
             rotate=180, vrange=None, colorbar=True, ax=None):
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
        rotate : int [180]
            Flag to rotate the mosaic by 90, 180, or 270 degrees.
            These rotations are implemented via a set of flips in x, y
            and/or xy-transpose operations.  If rotate is not in (0,
            90, 180, 270), a ValueError exception is raised.
        vrange : (float, float) [None]
            Range of pixel values to plot.  If None, then the cmap_range
            function will be used.
        colorbar : bool [True]
            Add a colorbar to the figure.
        ax : matplotlib.axes.Axes [None]
            Axes object to contain the figure. If None, then make a new
            figure.Figure object and add an axes.Axes object to use.
        """
        if ax is None:
            plt.rcParams['figure.figsize'] = figsize
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
        output_array = imutils.rebin_array(self.image_array, binsize,
                                           use_mean=True)
        if flipx:
            output_array = output_array[:, ::-1]
        if rotate not in (0, 90, 180, 270):
            raise ValueError(f'invalid rotation angle: {rotate}')
        if rotate == 180:
            ny, nx = output_array.shape
            rotated_array = np.zeros((nx, ny), dtype=output_array.dtype)
            for j in range(ny):
                rotated_array[:, ny-1-j] = output_array[::-1, j]
            output_array = rotated_array
        elif rotate == 90:
            rotated_array = copy.deepcopy(output_array.transpose())
            output_array = rotated_array[::-1, :]
        elif rotate == 270:
            rotated_array = copy.deepcopy(output_array.transpose())
            output_array = rotated_array[:, ::-1]
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
        if colorbar:
            fig.colorbar(image)
        # Turn off ticks and tick labels for x- and y-axes.
        plt.tick_params(axis='both', which='both',
                        top=False, bottom=False, left=False, right=False,
                        labelbottom=False, labelleft=False)
        # Label segments by sensor bay and segment number.
        for slot in self.fits_files:
            seg_coords = list(self._amp_coords[slot].values())[-8]
            xmin, xmax, ymin, ymax = seg_coords
            xx = float(xmax + xmin)/2./float(self.nx)
            if flipx:
                xx = 1 - xx
            yy = 1. - (float(ymax - ymin)*0.05 + ymin)/float(self.ny)
            if rotate == 180:
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
                if rotate == 180:
                    xx = 1 - xx
                    yy = 1 - yy
                plt.annotate('%s' % imutils.channelIds[amp],
                             (xx, yy), xycoords='axes fraction',
                             size='x-small', horizontalalignment='center',
                             verticalalignment='center', color=textcolor)
        plt.annotate(annotation, (1, -0.1), xycoords='axes fraction',
                     horizontalalignment='right', verticalalignment='bottom')


class CornerRaftMosaic:
    # The amp_llc data are the nominal locations in pixel coordinates
    # of the lower left corner of each imaging segment for the
    # corresponding amp in the guide and wavefront sensors.  These
    # values assume a coordinate system which is rotated so that the
    # wavefront sensors are in the lower left corner of the raft.
    # They are based on sheet 4 of LCA-13381.
    amp_llc = dict()
    amp_llc['SG0'] = {amp: (2050, 8226 - (amp - 1)*509) for amp in range(1, 9)}
    amp_llc['SG0'].update({amp: (49, 4663 + (amp - 9)*509)
                           for amp in range(9, 17)})
    amp_llc['SG1'] = {amp: (4663 + (amp - 1)*509, 2050) for amp in range(1, 9)}
    amp_llc['SG1'].update({amp: (8226 - (amp - 9)*509, 49)
                           for amp in range(9, 17)})
    amp_llc['SW0'] = {amp: (3777 - (amp-1)*509, 88) for amp in range(1, 9)}
    amp_llc['SW1'] = {amp: (214 + (amp-1)*509, 2413) for amp in range(1, 9)}
    wf_channels = {1: '10', 2: '11', 3: '12', 4: '13',
                   5: '14', 6: '15', 7: '16', 8: '17'}
    def __init__(self, fits_files, gains=None, bias_subtract=True,
                 nx=11630, ny=11630, bias_frames=None, dark_currents=None):
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
        nx : int [None]
            Number of pixels in the x (serial) direction.  If None,
            use 11630, which is based on the size of the corner raft
            baseplate in the x-direction.
        ny : int [None]
            Number of pixels in the y (parallel) direction.  If None,
            use 11630, which is based on the size of the corner raft
            baseplate in the y-direction.
        bias_frames : dict [None]
            Dictionary of single sensor bias frames, keyed by raft slot.
            If None, then just do the bias level subtraction using
            overscan region.
        dark_currents : dict [None]
            Dictionary of dictionaries of dark current values per amp
            in e-/s, keyed by raft slot and by amp number. If None, then
            dark current subtraction is not applied.
        """
        self.fits_files = fits_files
        with fits.open(list(fits_files.values())[0]) as hdus:
            self.raft_name = hdus[0].header['RAFTNAME']
            try:
                self.wl = hdus[0].header['MONOWL']
            except KeyError:
                self.wl = None
        if nx is None:
            nx = 11630
        if ny is None:
            ny = 11630
        self.image_array = np.zeros((nx, ny), dtype=np.float32)
        self.pixel_values = []
        self.nx = nx
        self.ny = ny
        if gains is None:
            # Assume unit gain for all amplifiers.
            unit_gains = dict([(i, 1) for i in range(1, 17)])
            gains = dict([(slot, unit_gains) for slot in fits_files])
        for slot, filename in fits_files.items():
            if slot not in self.amp_llc:
                continue
            bias_frame = bias_frames[slot] if bias_frames is not None else None
            ccd = sensorTest.MaskedCCD(filename, bias_frame=bias_frame)
            if dark_currents is not None:
                try:
                    dark_time = ccd.md.get('DARKTIME')
                except:
                    dark_time = ccd.md.get('EXPTIME')
            for amp in ccd:
                if dark_currents:
                    dark_correction = dark_time*dark_currents[slot][amp]
                else:
                    dark_correction = 0
                if slot.startswith('SW'):
                    self._set_sw_segment(slot, ccd, amp, gains[slot][amp],
                                         bias_subtract, dark_correction)
                else:
                    self._set_sg_segment(slot, ccd, amp, gains[slot][amp],
                                         bias_subtract, dark_correction)

    def _set_sw_segment(self, slot, ccd, amp, gain, bias_subtract,
                        dark_correction):
        if bias_subtract:
            mi = ccd.unbiased_and_trimmed_image(amp)
        else:
            mi = ccd[amp].Factory(ccd[amp], ccd.amp_geom.imaging)
        xmin, ymin = self.amp_llc[slot][amp]
        xmax = xmin + ccd.amp_geom.nx
        ymax = ymin + ccd.amp_geom.ny
        seg_array = np.array(gain*copy.deepcopy(mi.getImage().getArray()),
                             dtype=np.float32) - dark_correction
        self.pixel_values.extend(seg_array.ravel())
        if slot == 'SW0':
            # Flip amps in serial direction and parallel directions.
            # This is equivalent to a 180 degree rotation about the
            # segment center.
            seg_array = seg_array[::-1, ::-1]
        self.image_array[ymin:ymax, xmin:xmax] = seg_array

    def _set_sg_segment(self, slot, ccd, amp, gain, bias_subtract,
                        dark_correction):
        if bias_subtract:
            mi = ccd.unbiased_and_trimmed_image(amp)
        else:
            mi = ccd[amp].Factory(ccd[amp], ccd.amp_geom.imaging)
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

        seg_array = np.array(gain*copy.deepcopy(mi.getImage().getArray()),
                             dtype=np.float32) - dark_correction
        self.pixel_values.extend(seg_array.ravel())
        if amp < 9:
            seg_array = seg_array[::-1, :]
        if slot == 'SG0':
            seg_array = seg_array.transpose()
        else:
            seg_array = seg_array[:, ::-1]
        self.image_array[ymin:ymax, xmin:xmax] = seg_array

    def plot(self, title=None, cmap=plt.cm.hot, nsig=5, figsize=(10, 10),
             binsize=10, flipx=True, textcolor='c', annotation='',
             rotate=180, vrange=None, colorbar=True, ax=None):
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
        rotate : int [180]
            Flag to rotate the mosaic by 90, 180, or 270 degrees.
            These rotations are implemented via a set of flips in x, y
            and/or xy-transpose operations.  If rotate is not in (0,
            90, 180, 270), a ValueError exception is raised.
        vrange : (float, float) [None]
            Range of pixel values to plot.  If None, then the cmap_range
            function will be used.
        colorbar : bool [True]
            Add a colorbar to the figure.
        ax : matplotlib.axes.Axes [None]
            Axes object to contain the figure. If None, then make a new
            figure.Figure object and add an axes.Axes object to use.
        """
        if ax is None:
            plt.rcParams['figure.figsize'] = figsize
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
        output_array = imutils.rebin_array(self.image_array, binsize,
                                           use_mean=True)
        if flipx:
            output_array = output_array[:, ::-1]
        if rotate not in (0, 90, 180, 270):
            raise ValueError(f'invalid rotation angle: {rotate}')
        if rotate == 180:
            ny, nx = output_array.shape
            rotated_array = np.zeros((nx, ny), dtype=output_array.dtype)
            for j in range(ny):
                rotated_array[:, ny-1-j] = output_array[::-1, j]
            output_array = rotated_array
        elif rotate == 90:
            rotated_array = copy.deepcopy(output_array.transpose())
            output_array = rotated_array[::-1, :]
        elif rotate == 270:
            rotated_array = copy.deepcopy(output_array.transpose())
            output_array = rotated_array[:, ::-1]
        image = ax.imshow(output_array, interpolation='nearest', cmap=cmap)
        if vrange is None:
            # Set range and normalization of color map based on
            # sigma-clip of pixel values.
            num_pix = len(self.pixel_values)
            dpix = binsize*binsize
            rebinned_pixels = [sum(self.pixel_values[i*dpix:(i+1)*dpix])/dpix
                               for i in range(num_pix//dpix)]
            vmin, vmax = cmap_range(rebinned_pixels, nsig=nsig)
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
        if colorbar:
            fig.colorbar(image)
        # Turn off ticks and tick labels for x- and y-axes.
        plt.tick_params(axis='both', which='both',
                        top=False, bottom=False, left=False, right=False,
                        labelbottom=False, labelleft=False)
        if rotate != 180:
            # Only label segments if the corner raft orientation matches
            # sheet 4 of LCA-13381.
            return
        # Label segments by sensor bay and segment number.
        for slot in self.fits_files:
            if slot.startswith('SW'):
                channels = self.wf_channels
            else:
                channels = imutils.channelIds
            edge_offset = 350
            serial_midpoint = 509//2
            for amp, (x, y) in self.amp_llc[slot].items():
                if slot == 'SG0':
                    dx = edge_offset
                    dy = serial_midpoint
                    xx = (x + dx)/self.nx
                    if amp < 9:
                        xx += (2000 - 2*dx)/self.nx
                    yy = (y + dy)/self.ny
                else:
                    dx = serial_midpoint
                    dy = edge_offset
                    xx = (x + dx)/self.nx
                    yy = (y + dy)/self.ny
                    if ((slot == 'SG1' and amp < 9) or
                        slot == 'SW0'):
                        yy += (2000 - 2*dy)/self.ny
                plt.annotate(f'{channels[amp]}', (xx, yy),
                             xycoords='axes fraction',
                             horizontalalignment='center',
                             verticalalignment='center',
                             size='x-small', color=textcolor)
            if slot == 'SG1':
                x, y = self.amp_llc[slot][8]
                y += 2000 - edge_offset//2
            elif slot == 'SG0':
                x, y = self.amp_llc[slot][9]
                y += (509 - edge_offset//2)
            elif slot == 'SW1':
                x, y = self.amp_llc[slot][1]
                y += 2000 - edge_offset//2
            elif slot == 'SW0':
                x, y = self.amp_llc[slot][8]
                y += 2000 - edge_offset//2
            plt.annotate(f'{slot}', (x/self.nx, y/self.ny),
                         xycoords='axes fraction', size='x-small',
                         horizontalalignment='left',
                         verticalalignment='bottom',
                         color=textcolor)
