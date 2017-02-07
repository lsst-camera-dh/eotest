"""
Code to perform raft-level mosaicking from single sensor frames
compliant with LCA-13501.
"""
from __future__ import absolute_import, print_function
import os
import copy
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
    def __init__(self, fits_files, gains_list=None, bias_subtract=True,
                 nx=12700, ny=12700):
        """
        Constructor.

        Parameters
        ----------
        fits_files : list
            Nine item list of single sensor FITS files with header
            keywords conforming to LCA-13501.
        gains_list : list, optional
            Nine item list of dictionaries (one per FITS file) of
            system gain values for each amp.  Default: None (do not
            apply gain correction).
        bias_subtract : bool, optional
            Flag do to a bias subtraction based on the serial overscan.
            Default: True
        nx : int, optional
            Number of pixels in the x (serial) direction.  Default: 12700.
        ny : int, optional
            Number of pixels in the y (parallel) direction.  Default: 12700.
        """
        self.raft_name = fits.open(fits_files[0])[0].header['RAFTNAME']
        self.image_array = np.zeros((nx, ny), dtype=np.float32)
        if gains_list is None:
            # Assume unit gain for all amplifiers.
            unit_gains = dict([(i, 1) for i in range(1, 17)])
            gains_list = [unit_gains]*len(fits_files)

        for item, ccd_gains in zip(fits_files, gains_list):
            print("processing", os.path.basename(item))
            ccd = sensorTest.MaskedCCD(item)
            hdu_list = fits.open(item)
            for amp, hdu in zip(ccd, hdu_list[1:]):
                amp_gain = ccd_gains[amp]
                self._set_segment(ccd, amp, hdu, amp_gain, bias_subtract)

    def _set_segment(self, ccd, amp, hdu, amp_gain, bias_subtract):
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
        # Write the segment pixel values into the full raft image mosaic.
        self.image_array[ymin:ymax, xmin:xmax] = seg_array

    def plot(self, cmap=plt.cm.hot, nsig=5, figsize=(10, 10)):
        """
        Render the raft mosaic.
        """
        plt.rcParams['figure.figsize'] = figsize
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        image = ax.imshow(self.image_array, interpolation='nearest', cmap=cmap)
        vmin, vmax = cmap_range(self.image_array, nsig=nsig)
        norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
        image.set_norm(norm)
        ax.set_title(self.raft_name)
        fig.colorbar(image)
        return fig
