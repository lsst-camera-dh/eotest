"""
Code to perform raft-level mosaicking from single sensor frames
compliant with LCA-13501.
"""
from __future__ import absolute_import, print_function
import os
import copy
from collections import defaultdict
import numpy as np
import astropy.io.fits as fits
import matplotlib
import matplotlib.pyplot as plt
import lsst.eotest.sensor as sensorTest
from lsst.eotest.sensor.EOTestPlots import cmap_range

__all__ = ['RaftMosaic']

class RaftMosaic(object):
    """
    Raft level mosaic of individual CCDs.
    """
    _segment_mapping = dict([(1, '10'), (2, '11'), (3, '12'), (4, '13'),
                             (5, '14'), (6, '15'), (7, '16'), (8, '17'),
                             (9, '07'), (10, '06'), (11, '05'), (12, '04'),
                             (13, '03'), (14, '02'), (15, '01'), (16, '00')])
    def __init__(self, fits_files, gains=None, bias_subtract=True,
                 nx=12700, ny=12700, nx_segments=8, ny_segments=2):
        """
        Constructor.

        Parameters
        ----------
        fits_files : dict
            Dictionary of single sensor FITS files, keyed by raft slot
            name.  These files should conform to LCA-13501.
        gains : dict, optional
            Dictionary (keyed by slot name) of dictionaries (one per
            FITS file) of system gain values for each amp.  Default:
            None (i.e., do not apply gain correction).
        bias_subtract : bool, optional
            Flag do to a bias subtraction based on the serial overscan.
            Default: True
        nx : int, optional
            Number of pixels in the x (serial) direction.  Default: 12700
        ny : int, optional
            Number of pixels in the y (parallel) direction.  Default: 12700
        nx_segments : int, optional
            Number of segments in the x (serial) direction.  Default: 8
        ny_segments : int, optional
            Number of pixels in the y (parallel) direction.  Default: 2
        """
        self.fits_files = fits_files
        self.raft_name = fits.open(fits_files.values()[0])[0].header['RAFTNAME']
        self.wl = fits.open(fits_files.values()[0])[0].header['MONOWL']
        self.image_array = np.zeros((nx, ny), dtype=np.float32)
        self.nx = nx
        self.ny = ny
        self.nx_segments = nx_segments
        self.ny_segments = ny_segments
        self._amp_coords = defaultdict(dict)
        if gains is None:
            # Assume unit gain for all amplifiers.
            unit_gains = dict([(i, 1) for i in range(1, 17)])
            gains = dict([(slot, unit_gains) for slot in fits_files])
        for slot, filename in fits_files.items():
            print("processing", os.path.basename(filename))
            ccd = sensorTest.MaskedCCD(filename)
            hdu_list = fits.open(filename)
            for amp, hdu in zip(ccd, hdu_list[1:]):
                self._set_segment(slot, ccd, amp, hdu, gains[slot][amp],
                                  bias_subtract)

    def _set_segment(self, slot, ccd, amp, hdu, amp_gain, bias_subtract):
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
        # Determine flip in serial direction based on 1, 1 element of
        # transformation matrix.
        if hdu.header['PC1_1Q'] < 0:
            seg_array = seg_array[::-1, :]
            xmax = int(hdu.header['CRVAL1Q'])
            xmin = xmax - ccd.amp_geom.nx
        else:
            xmin = int(hdu.header['CRVAL1Q'])
            xmax = xmin + ccd.amp_geom.nx
        # Determine flip in parallel direction based on 2, 2 element
        # of transformation matrix.
        if hdu.header['PC2_2Q'] < 0:
            seg_array = seg_array[:, ::-1]
            ymax = int(hdu.header['CRVAL2Q'])
            ymin = ymax - ccd.amp_geom.ny
        else:
            ymin = int(hdu.header['CRVAL2Q'])
            ymax = ymin + ccd.amp_geom.ny
        # Save coordinates of segment for later use.
        self._amp_coords[slot][amp] = xmin, xmax, ymin, ymax
        # Write the segment pixel values into the full raft image mosaic.
        self.image_array[ymin:ymax, xmin:xmax] = seg_array

    def plot(self, title=None, cmap=plt.cm.hot, nsig=5, figsize=(10, 10)):
        """
        Render the raft mosaic.
        """
        plt.rcParams['figure.figsize'] = figsize
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        image = ax.imshow(self.image_array, interpolation='nearest', cmap=cmap)
        # Set range and normalization of color map based on sigma-clip
        # of pixel values.
        vmin, vmax = cmap_range(self.image_array, nsig=nsig)
        norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
        image.set_norm(norm)
        if title is None:
            title = "%s, %i nm" % (self.raft_name, self.wl)
        ax.set_title(title)
        fig.colorbar(image)
        # Turn off tick labels for x- and y-axes.
        plt.setp(ax.get_xticklabels(), visible=False)
        plt.setp(ax.get_yticklabels(), visible=False)
        # Label segments by sensor bay and segment number.
        for slot in self.fits_files:
            seg_coords = self._amp_coords[slot].values()[0]
            xmin, xmax, ymin, ymax = seg_coords
            xx = float(xmax + xmin)/2./float(self.nx)
            yy = (float(ymax - ymin)*0.95 + ymin)/float(self.ny)
            plt.annotate('%s' % slot,
                         (xx, yy), xycoords='axes fraction',
                         size='x-small', horizontalalignment='center',
                         verticalalignment='center')
            for amp, seg_coords in self._amp_coords[slot].items():
                xmin, xmax, ymin, ymax = seg_coords
                xx = float(xmax + xmin)/2./float(self.nx)
                if amp <= 8:
                    yy = (float(ymax - ymin)*0.85 + ymin)/float(self.ny)
                if amp > 8:
                    yy = (float(ymax - ymin)*0.15 + ymin)/float(self.ny)
                plt.annotate('%s' % self._segment_mapping[amp],
                             (xx, yy), xycoords='axes fraction',
                             size='x-small', horizontalalignment='center',
                             verticalalignment='center')
        return fig
