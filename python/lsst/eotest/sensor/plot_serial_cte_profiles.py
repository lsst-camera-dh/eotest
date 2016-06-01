"""
Plot the profile of mean overscan columns as a function of column # to
illustrate the serial CTE.
"""
from __future__ import print_function, absolute_import, division
import numpy as np
import matplotlib.pyplot as plt
import lsst.afw.math as afwMath
import lsst.eotest.sensor as sensorTest
from lsst.eotest.Estimator import Estimator

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

def bias_estimate(masked_image, amp_geom, nskip_last_cols=4):
    """
    Estimate the bias level and error in serial overscan region,
    skipping the nskip_last_cols columns to avoid bright columns in
    the e2v data.
    """
    imarr = masked_image.getImage().getArray()
    overscan = imarr[:amp_geom.serial_overscan.getMaxY()+1,
                     amp_geom.serial_overscan.getMinX()+overscans:-nskip_last_cols]
    bias_est = Estimator()
    bias_stats = afwMath.makeStatistics(np.array(overscan.flatten(),
                                                 dtype=np.float),
                                        afwMath.MEAN | afwMath.STDEV)
    bias_est.value = bias_stats.getValue(afwMath.MEAN)
    bias_est.error = bias_stats.getValue(afwMath.STDEV)
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

def get_overscan_columns(masked_image, gain, amp_geom, bias_est=None):
    "Compute the column Estimators, including bias subtraction"
    if bias_est is None:
        bias_est = gain*bias_estimate(masked_image, amp_geom)
    imarr = gain*masked_image.getImage().getArray()
    num_columns = amp_geom.full_segment.getMaxX()
    ymax = amp_geom.imaging.getMaxY()
    estimators = EstimatorList([make_estimator(imarr[:ymax+1, ix-1], bias_est)
                                for ix in range(1, num_columns+1)])
    return estimators

def plot_cte_profile(axes, masked_image, gain, amp_geom, cti, bias_est=None,
                     bias_subtract=True):
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
    column_ests = get_overscan_columns(masked_image, gain, amp_geom,
                                       bias_est=bias_est)
    columns = np.arange(1, amp_geom.full_segment.getMaxX()+1, dtype=np.float)
    plt.step(columns, column_ests.values-bias_offset, where='mid', color=color)
    axisrange = list(plt.axis())
    axisrange[:2] = (amp_geom.imaging.getMaxX() + 2 - 0.5,
                     amp_geom.full_segment.getMaxX() + 0.5)
    axisrange[2:] = (min(column_ests.values[amp_geom.imaging.getMaxX() + 1:])
                     - bias_offset - 3,
                     max(column_ests.values[amp_geom.imaging.getMaxX() + 1:])
                     - bias_offset + 3)
    plt.axis(axisrange)
    axes.errorbar(columns, column_ests.values-bias_offset,
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

    axes.errorbar([axisrange[0] + float(overscans)/2.], [pred_trail.value],
                  xerr=float(overscans)/2., yerr=pred_trail.error, fmt='o',
                  color=pred_color)

def plot_serial_cte_profiles(ccd, gains, cti, bias_est=None, figsize=(10, 10)):
    """
    Plot a 4x4 array of serial cte profiles for all 16 channels.
    """
    if bias_est is None:
        bias_est = dict((amp, gains[amp]*bias_estimate(ccd[amp], ccd.amp_geom))
                        for amp in ccd)
    plt.rcParams['figure.figsize'] = figsize
    fig = plt.figure()
    frame_axes = fig.add_subplot(111, frame_on=False, xticklabels=(),
                                 yticklabels=())
    frame_axes.set_xlabel('\ncolumn #')
    frame_axes.set_ylabel('column mean - bias (e-/pixel)\n')
    for amp in ccd:
        subplot = (4, 4, amp)
        axes = fig.add_subplot(*subplot)
        plot_cte_profile(axes, ccd[amp], gains[amp], ccd.amp_geom,
                         cti[amp], bias_est[amp], bias_subtract=True)
    return fig, frame_axes

if __name__ == '__main__':
    overscans = 2
    sensor_id = 'ITL-3800C-007'
#    flux_level = 'low'
    flux_level = 'high'
    sflat = '%s_superflat_%s.fits' % (sensor_id, flux_level)
    results = sensorTest.EOTestResults('%s_eotest_results.fits' % sensor_id)
    gains = dict((amp, gain) for amp, gain
                 in zip(results['AMP'], results['GAIN']))
    ccd = sensorTest.MaskedCCD(sflat)
    cti = dict([(amp, Estimator()) for amp in ccd])
    for amp in cti:
        cti[amp].value = results['CTI_%s_SERIAL' % flux_level.upper()][amp-1]
        cti[amp].error = results['CTI_%s_SERIAL_ERROR' % flux_level.upper()][amp-1]

    fig, frame_axes = plot_serial_cte_profiles(ccd, gains, cti)
    frame_axes.set_title(sflat)
