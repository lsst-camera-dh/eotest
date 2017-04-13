"""
A function to make histograms of total noise, particularly including
the dark current contribution, for specified exposure times.
"""
from __future__ import absolute_import
__all__ == ['total_noise_histogram']
import numpy as np
import matplotlib.pyplot as plt
import lsst.eotest.image_utils as imutils

def total_noise_histograms(dark_curr_pixels, read_noise, dark95s, exptime=16,
                           title=None, figsize=(10, 10)):
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
    plt.pcParams['figure.figsize'] = figsize
    fig = plt.figure()
    frame_axes = fig.add_subplot(111, frameon=False)
    if title is not None:
        frame_axes.set_title(title)
    frame_axes.set_xlabel('total noise (rms e-/pixel)')
    frame_axes.get_xaxis().set_ticks([])
    frame_axes.get_yaxis().set_ticks([])
    for amp in dark_curr_pixels:
        subplot = (4, 4, amp)
        ax = fig.add_subplot(*subplot)
        noise_values \
            = np.sqrt(dark_curr_pixels[amp]*exptime + read_noise[amp]**2)
        noise_dc95 = np.sqrt(dark95s[amp]*exptime + read_noise[amp]**2)
        plt.hist(noise_values, bins=100, histtype='step', color='blue')
        xmin, xmax, ymin, ymax = plt.axis()
        plt.plot([noise_dc95, noise_dc95], [ymin, ymax], ':k',
                 markersize=3, linewidth=1)
        plt.annotate('Segment %s' % imutils.channelIds[amp], (0.5, 0.9),
                     xycoords='axes fraction', size='x-small')
    return fig
