from __future__ import print_function
import os
import numpy as np
import matplotlib.pyplot as plt
import lsst.afw.math as afwMath
import lsst.eotest.sensor as sensorTest
from lsst.eotest.Estimator import Estimator

plt.ion()

def make_estimator(y, bias):
    estimator = Estimator()
    yvals = np.array(y.flatten(), dtype=np.float)
    estimator.value = np.mean(yvals)
    poisson_variance_per_pixel = np.sum(yvals - bias.value)/len(yvals)**2
    if poisson_variance_per_pixel > 0:
        estimator.error = np.sqrt(poisson_variance_per_pixel + bias.error**2)
    else:
        estimator.error = bias.error
    return estimator

def bias_estimate(masked_image, amp_geom):
    imarr = masked_image.getImage().getArray()
    overscan = imarr[:amp_geom.serial_overscan.getMaxY()+1,
                      amp_geom.serial_overscan.getMinX()+overscans:-4]
    bias_est = Estimator()
    bias_stats = afwMath.makeStatistics(np.array(overscan.flatten(),
                                                 dtype=np.float),
                                        afwMath.MEANCLIP | afwMath.STDEV)
    bias_est.value = bias_stats.getValue(afwMath.MEANCLIP)
    bias_est.error = bias_stats.getValue(afwMath.STDEV)
    return bias_est

class EstimatorList(list):
    "Container list for Estimator objects"
    def __init__(self, *args, **kwds):
        "Constructor"
        super(EstimatorList, self).__init__(*args, **kwds)

    @property
    def values(self):
        "Return a list of the values."
        return [x.value for x in self]

    @property
    def errors(self):
        "Return a list of the errors."
        return [x.error for x in self]

def get_overscan_columns(masked_image, gain, amp_geom, bias_est=None):
    "Compute the column Estimators, including bias subtraction"
    if bias_est is None:
        bias_est = gain*bias_estimate(masked_image, amp_geom)
    imarr = gain*masked_image.getImage().getArray()
    columns = range(1, amp_geom.full_segment.getMaxX() + 1)
    ymax = amp_geom.imaging.getMaxY()
    estimators = EstimatorList([make_estimator(imarr[:ymax+1, ix-1], bias_est)
                                for ix in columns])
    return estimators

def plot_cte_profile(axes, masked_image, gain, amp_geom, cti, bias_est,
                     bias_subtract=True):
    if bias_subtract:
        bias_offset = bias_est.value
    else:
        bias_offset = 0

    color = 'blue'
    pred_color = 'red'
    column_ests = get_overscan_columns(masked_image, gain, amp_geom,
                                       bias_est=bias_est)
    columns = range(1, amp_geom.full_segment.getMaxX() + 1)
    plt.step(columns, column_ests.values-bias_offset, where='mid', color=color)
    axisrange = list(plt.axis())
    axisrange[:2] = (amp_geom.imaging.getMaxX() + 2-0.5,
                     amp_geom.imaging.getMaxX() + 24)
    axisrange[2:] = (min(column_ests.values[amp_geom.imaging.getMaxX() + 1:])
                     - bias_offset - 3,
                     max(column_ests.values[amp_geom.imaging.getMaxX() + 1:])
                     - bias_offset + 3)
    plt.axis(axisrange)
    axes.errorbar(columns, column_ests.values-bias_offset,
                  yerr=column_ests.errors,
                  xerr=0.5, fmt='o', color=color)

    ntransfers = amp_geom.imaging.getMaxX()
    Si = column_ests[amp_geom.imaging.getMaxX()]

    if bias_subtract:
        plt.plot(axisrange[:2], (0, 0), 'k:')
    else:
        plt.plot(axisrange[:2], (bias_est.value, bias_est.value), 'k:')

    pred_trail = (cti*ntransfers*(Si - bias_est.value)/float(overscans)
                  + bias_est) - bias_offset

    axes.errorbar([axisrange[0] + overscans/2.], [pred_trail.value],
                  xerr=overscans/2., yerr=pred_trail.error, fmt='o',
                  color=pred_color)

def plot_serial_cte_profiles(ccd, gains, cti, bias_est, figsize=(10, 10)):
    plt.rcParams['figure.figsize'] = figsize
    fig  = plt.figure()
    frame_axes = fig.add_subplot(111, frame_on=False, xticklabels=(),
                                 yticklabels=())
    frame_axes.set_xlabel('\ncolumn #')
    frame_axes.set_ylabel('Mean e-/pixel - bias\n')
    for amp in ccd:
        subplot = (4, 4, amp)
        axes = fig.add_subplot(*subplot)
        plot_cte_profile(axes, ccd[amp], gains[amp], ccd.amp_geom,
                         cti[amp], bias_est[amp], bias_subtract=True)
    return fig, frame_axes

if __name__ == '__main__':
    overscans = 2
    sflat = 'E2V-CCD250-114_superflat_low.fits'
    #sflat = 'E2V-CCD250-114_superflat_high.fits'
    results = sensorTest.EOTestResults('E2V-CCD250-114_eotest_results.fits')
    gains = dict((amp, gain) for amp, gain
                 in zip(results['AMP'], results['GAIN']))
    ccd = sensorTest.MaskedCCD(sflat)
    #
    # run EPERTask
    #
    eper_task = sensorTest.EPERTask()
    eper_task.log.setThreshold(eper_task.log.INFO)
    eper_task.config.verbose = False
    eper_task.config.cti = True
    eper_task.config.direction = 's'
    cti, bias_est = eper_task.run(sflat, range(1, 17), overscans, gains=gains)
    for amp in cti:
        print(amp, cti[amp])

    fig, frame_axes = plot_serial_cte_profiles(ccd, gains, cti, bias_est)
    frame_axes.set_title(sflat)
