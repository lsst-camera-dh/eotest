"""
A function to make histograms of total noise, particularly including
the dark current contribution, for specified exposure times.
"""
from __future__ import absolute_import
import numpy as np
import matplotlib.pyplot as plt
import lsst.afw.math as afwMath
import lsst.eotest.image_utils as imutils

__all__ = ['total_noise_histograms']

def total_noise_histograms(dark_curr_pixels, read_noise, dark95s, exptime=16,
                           title=None, bins=50, figsize=(10, 10)):
    """
    Make histograms of total noise per pixel. Total noise is defined
    as the device read and system noise for each amp combined in
    quadrature with the inferred Poisson noise from dark current for a
    given exposure time.

    Parameters
    ----------
    dark_curr_pixels : dict
        A dictionary of np arrays of dark current values for unmasked
        pixels, keyed by amp number.
    read_noise : dict
        A dictionary of device combined read and system noise, keyed by
        amp number. Units are rms e-/pixel.
    dark95s : dict
        The 95th percentile dark current pixel value for each amp, keyed
        by amp.  Units are e-/s.
    exptime : float, optional
        The exposure time in seconds to use for computing the dark current
        contribution to the per pixel noise.  Default: 16.
    title : str, optional
        The title of the figure. If None (default), do not set a title.
    figsize : (float, float), optional
        The overall size (in inches) of the 4x4 array of histograms.
        Default: (10, 10)

    Returns
    -------
    matplotlib.figure.Figure
        The figure with the histograms of total noise for each amp.
    """
    plt.rcParams['figure.figsize'] = figsize
    fig = plt.figure()
    frame_axes = fig.add_subplot(111, frameon=False)
    if title is not None:
        frame_axes.set_title(title)
    frame_axes.set_xlabel('\ntotal noise (rms e-/pixel)')
    frame_axes.set_ylabel('normalized histogram\n')
    frame_axes.get_xaxis().set_ticks([])
    frame_axes.get_yaxis().set_ticks([])
    for amp in dark_curr_pixels:
        subplot = fig.add_subplot(4, 4, amp)
        subplot.tick_params(axis='both', labelsize='x-small')
        noise_values \
            = np.sqrt(dark_curr_pixels[amp]*exptime + read_noise[amp]**2)
        stats = afwMath.makeStatistics(np.array(noise_values, dtype=np.float),
                                       afwMath.MEDIAN | afwMath.STDEVCLIP)
        median = stats.getValue(afwMath.MEDIAN)
        stdev = stats.getValue(afwMath.STDEVCLIP)
        x_range = (median - 5.*stdev, median + 5.*stdev)
        noise_dc95 = np.sqrt(dark95s[amp]*exptime + read_noise[amp]**2)
        plt.hist(noise_values, bins=bins, histtype='step', color='blue',
                 range=x_range, normed=True)
        ymin, ymax = plt.axis()[2:]
        plt.plot([noise_dc95, noise_dc95], [ymin, ymax], ':k',
                 markersize=3, linewidth=1)
        plt.annotate('Segment %s' % imutils.channelIds[amp], (0.5, 0.9),
                     xycoords='axes fraction', size='x-small')
    plt.tight_layout()
    return fig, frame_axes
