"""
Plot the profile of mean overscan columns as a function of column # to
illustrate the serial CTE.
"""
from __future__ import print_function, absolute_import, division
import numpy as np
import matplotlib.pyplot as plt
from . import pylab_plotter as plot
import lsst.afw.math as afwMath
from lsst.eotest.Estimator import Estimator

__all__ = ['bias_estimate', 'serial_cte_profile', 'parallel_cte_profile']

plt.ion()

def make_estimator(pixel_values, bias):
    """
    Create an Estimator object (value and error) based on the counts
    in the imaging region column after bias subtraction. Gain is
    assumed to have been applied to both the column pixel values
    and the bias level.
    """
    estimator = Estimator()
    yvals = np.array(pixel_values.flatten(), dtype=np.float)
    estimator.value = np.mean(yvals)
    poisson_variance_per_pixel = np.sum(yvals - bias.value)/float(len(yvals)**2)
    if poisson_variance_per_pixel > 0:
        estimator.error = np.sqrt(poisson_variance_per_pixel + bias.error**2)
    else:
        estimator.error = bias.error
    return estimator

def bias_estimate(masked_image, amp_geom, overscans=2, nskip_last_cols=4,
                  serial=True):
    """
    Estimate the bias level and error in serial or parallel overscan
    region, skipping the nskip_last_cols columns to avoid bright
    columns in the e2v data.
    """
    imarr = masked_image.getImage().getArray()
    if serial:
        overscan = imarr[:amp_geom.serial_overscan.getMaxY() + 1,
                          amp_geom.serial_overscan.getMinX() + overscans:
                              -nskip_last_cols].flatten()
    else:
        overscan = imarr[amp_geom.parallel_overscan.getMinY() + overscans:,
                         :amp_geom.parallel_overscan.getMaxX() + 1].flatten()
    bias_est = Estimator()
    bias_stats = afwMath.makeStatistics(np.array(overscan, dtype=np.float),
                                        afwMath.MEAN | afwMath.STDEV)
    bias_est.value = bias_stats.getValue(afwMath.MEAN)
    bias_est.error = bias_stats.getValue(afwMath.STDEV)/np.sqrt(len(overscan))
    return bias_est

class EstimatorList(list):
    "Container list for Estimator objects"
    def __init__(self, *args, **kwds):
        "Constructor"
        super(EstimatorList, self).__init__(*args, **kwds)

    @property
    def values(self):
        "Return a numpy.array of the values."
        return np.array([x.value for x in self])

    @property
    def errors(self):
        "Return a numpy.array of the errors."
        return np.array([x.error for x in self])

def get_overscan_ests(masked_image, gain, amp_geom, bias_est=None,
                      overscans=2, serial=True):
    "Compute the column Estimators, including bias subtraction"
    if bias_est is None:
        bias_est = gain*bias_estimate(masked_image, amp_geom,
                                      overscans=overscans, serial=serial)
    imarr = gain*masked_image.getImage().getArray()
    if serial:
        num_columns = amp_geom.full_segment.getMaxX()
        ymax = amp_geom.imaging.getMaxY()
        estimators = EstimatorList([make_estimator(imarr[:ymax+1, ix-1],
                                                   bias_est)
                                    for ix in range(1, num_columns+1)])
    else:
        num_rows = amp_geom.full_segment.getMaxY()
        xmax = amp_geom.imaging.getMaxX()
        estimators = EstimatorList([make_estimator(imarr[iy-1, :xmax+1],
                                                   bias_est)
                                    for iy in range(1, num_rows+1)])
    return estimators

def serial_cte_profile(axes, masked_image, gain, amp_geom, cti, bias_est,
                       bias_subtract=True, overscans=2, xaxis_range=None,
                       cti_spec=5e-6):
    """
    Plot the gain-corrected mean column pixel values vs column number.
    Overlay the mean trailed signal in the selected overscan columns.
    """
    color = 'blue'
    pred_color = 'red'
    if bias_subtract:
        bias_offset = bias_est.value
    else:
        bias_offset = 0
    column_ests = get_overscan_ests(masked_image, gain, amp_geom,
                                    bias_est=bias_est, overscans=overscans)
    columns = np.arange(1, amp_geom.full_segment.getMaxX()+1, dtype=np.float)
    plt.step(columns, column_ests.values-bias_offset, where='mid', color=color)
    axisrange = list(plt.axis())
    if xaxis_range is None:
        axisrange[:2] = (amp_geom.imaging.getMaxX() + 2 - 0.5,
                         amp_geom.full_segment.getMaxX() + 0.5)
        axisrange[2:] = (min(column_ests.values[amp_geom.imaging.getMaxX()+1:])
                         - bias_offset - 10,
                         max(column_ests.values[amp_geom.imaging.getMaxX()+1:])
                         - bias_offset + 10)
    else:
        axisrange[:2] = xaxis_range
    plt.axis(axisrange)
    axes.errorbar(columns, column_ests.values - bias_offset,
                  yerr=column_ests.errors,
                  xerr=0.5, fmt='o', color=color)

    ntransfers = amp_geom.imaging.getMaxX()
    si = column_ests[amp_geom.imaging.getMaxX()]

    if bias_subtract:
        plt.plot(axisrange[:2], (0, 0), 'k:')
    else:
        plt.plot(axisrange[:2], (bias_est.value, bias_est.value), 'k:')

    pred_trail = (cti*ntransfers*(si - bias_est.value)/float(overscans)
                  + bias_est) - bias_offset
    axes.errorbar([amp_geom.imaging.getMaxX() + 1.5 + float(overscans)/2.],
                  [pred_trail.value],
                  xerr=float(overscans)/2., yerr=pred_trail.error, fmt='o',
                  color=pred_color)
    spec = (cti_spec*ntransfers*(si - bias_est.value)/float(overscans)
            + bias_est) - bias_offset
    plt.plot(axisrange[:2], (spec.value, spec.value), 'r:')

def parallel_cte_profile(axes, masked_image, gain, amp_geom, cti, bias_est,
                         bias_subtract=True, overscans=2, xaxis_range=None,
                         cti_spec=3e-6):
    """
    Plot the gain-corrected mean row pixel values vs row number.
    Overlay the mean trailed signal in the selected overscan rows.
    """
    color = 'blue'
    pred_color = 'red'
    if bias_subtract:
        bias_offset = bias_est.value
    else:
        bias_offset = 0
    row_ests = get_overscan_ests(masked_image, gain, amp_geom,
                                 bias_est=bias_est, overscans=overscans,
                                 serial=False)
    rows = np.arange(1, amp_geom.full_segment.getMaxY()+1, dtype=np.float)
    plt.step(rows, row_ests.values-bias_offset, where='mid', color=color)
    axisrange = list(plt.axis())
    if xaxis_range is None:
        axisrange[:2] = (amp_geom.imaging.getMaxY() + 2 - 0.5,
                         amp_geom.full_segment.getMaxY() + 0.5)
        axisrange[2:] = (min(row_ests.values[amp_geom.imaging.getMaxY()+1:])
                         - bias_offset - 10,
                         max(row_ests.values[amp_geom.imaging.getMaxY()+1:])
                         - bias_offset + 10)
    else:
        axisrange[:2] = xaxis_range
    plt.axis(axisrange)
    axes.errorbar(rows, row_ests.values - bias_offset,
                  yerr=row_ests.errors, xerr=0.5, fmt='o', color=color)

    ntransfers = amp_geom.imaging.getMaxY()
    si = row_ests[amp_geom.imaging.getMaxY()]

    if bias_subtract:
        plt.plot(axisrange[:2], (0, 0), 'k:')
    else:
        plt.plot(axisrange[:2], (bias_est.value, bias_est.value), 'k:')

    pred_trail = (cti*ntransfers*(si - bias_est.value)/float(overscans)
                  + bias_est) - bias_offset
    axes.errorbar([amp_geom.imaging.getMaxY() + 1.5 + float(overscans)/2.],
                  [pred_trail.value],
                  xerr=float(overscans)/2., yerr=pred_trail.error, fmt='o',
                  color=pred_color)
    spec = (cti_spec*ntransfers*(si - bias_est.value)/float(overscans)
            + bias_est) - bias_offset
    plt.plot(axisrange[:2], (spec.value, spec.value), 'r:')

def plot_cte_profiles(ccd, gains, cti, title=None, bias_est=None,
                      bias_subtract=True,
                      xaxis_range=None, figsize=(11, 8.5), serial=True):
    """
    Plot a 4x4 array of serial or parallel cte profiles for all 16 channels.
    """
    if bias_est is None:
        bias_est = dict((amp, gains[amp]*bias_estimate(ccd[amp], ccd.amp_geom,
                                                       serial=serial))
                        for amp in ccd)
    if serial:
        xlabel = 'column #'
    else:
        xlabel = 'row #'
    if bias_subtract:
        ylabel = 'mean signal - bias (e-/pixel)'
    else:
        ylabel = 'mean signal (e-/pixel)'
    xoffset, yoffset = 0.025, 0.025
    for amp in ccd:
        subplot = (4, 4, amp)
        if amp == 1:
            win = plot.Window(subplot=subplot, figsize=figsize,
                              xlabel=xlabel, ylabel=ylabel, size='large')
            if title is not None:
                win.frameAxes.text(0.5, 1.08, title,
                                   horizontalalignment='center',
                                   verticalalignment='top',
                                   transform=win.frameAxes.transAxes,
                                   size='large')
        else:
            win.select_subplot(*subplot)
        bbox = win.axes[-1].get_position()
        points = bbox.get_points()
        points[0] += xoffset
        points[1] += yoffset
        bbox.set_points(points)
        win.axes[-1].set_position(bbox)
        if serial:
            serial_cte_profile(win.axes[-1], ccd[amp], gains[amp],
                               ccd.amp_geom, cti[amp], bias_est[amp],
                               bias_subtract=bias_subtract,
                               xaxis_range=xaxis_range)
        else:
            parallel_cte_profile(win.axes[-1], ccd[amp], gains[amp],
                                 ccd.amp_geom, cti[amp], bias_est[amp],
                                 bias_subtract=bias_subtract,
                                 xaxis_range=xaxis_range)
    return win

if __name__ == '__main__':
    from .EOTestResults import EOTestResults
    from .MaskedCCD import MaskedCCD
    sensor_id = 'ITL-3800C-012'
#    flux_level = 'low'
    flux_level = 'high'
    sflat = '%s_superflat_%s.fits' % (sensor_id, flux_level)
    results = EOTestResults('%s_eotest_results.fits' % sensor_id)
    gains = dict((amp, gain) for amp, gain
                 in zip(results['AMP'], results['GAIN']))
    ccd = MaskedCCD(sflat)
    cti = dict([(amp, Estimator()) for amp in ccd])
    for amp in cti:
        cti[amp].value = results['CTI_%s_SERIAL' % flux_level.upper()][amp-1]
        cti[amp].error = results['CTI_%s_SERIAL_ERROR' % flux_level.upper()][amp-1]
    fig, frame_axes = plot_serial_cte_profiles(ccd, gains, cti)
    frame_axes.set_title(sflat)
