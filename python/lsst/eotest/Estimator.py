from __future__ import division
import numpy as np
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath

class Estimator(object):
    "Abstraction for a point estimator of pixel data and its errors"
    def __init__(self, *args, **kwds):
        self.image = None
        self.gain = None
        self.stat_ctrl = None
        self.statistic = None
        self.value = None
        self.error = None
        self._format_str = None
        if args:
            self.set_properties(*args, **kwds)
    def set_properties(self, image, stat_ctrl, gain=1, statistic=afwMath.MEAN,
                       var_wt=1):
        # Make a deep copy of the input image so that we can convert to
        # e- and have Poisson statistics apply.
        self.image = image.clone()
        self.image *= gain
        self.stat_ctrl = stat_ctrl
        self.gain = gain
        self.statistic = statistic
        self.var_wt = var_wt
        self._compute_stats()
    def _compute_stats(self):
        if self.stat_ctrl is None:
            makeStatistics = lambda *args : afwMath.makeStatistics(*args[:2])
        else:
            makeStatistics = lambda *args : afwMath.makeStatistics(*args[:3])
        if self.statistic not in (afwMath.MEAN, afwMath.MEDIAN):
            # In case other statistics are given, set error to zero for now.
            self.value = makeStatistics(self.image, self.statistic,
                                        self.stat_ctrl).getValue()
            self.error = 0
            return
        # Compute the error assuming the statistic is afw.MEAN.  For
        # Gaussian stats, the error on the median is sqrt(pi/2.)
        # times the error on the mean, but for Poisson stats, it is
        # actually zero when the number of pixels is much larger than
        # the expected count per pixel, but within factors of order
        # unity to the error on the mean for numpix \la O(100)*count/pixel.
        flags = self.statistic | afwMath.SUM | afwMath.MEAN
        stats = makeStatistics(self.image, flags, self.stat_ctrl)
        pixel_sum = stats.getValue(afwMath.SUM)
        # Infer the number of pixels taking into account masking.
        if pixel_sum == 0:
            # Handle case where pixel_sum is zero (and hence the
            # mean is zero).
            self.value = 0
            self.error = 0
            return
        npix = pixel_sum/stats.getValue(afwMath.MEAN)
        self.value = stats.getValue(self.statistic)
        self.error = np.sqrt(pixel_sum/npix/self.var_wt)
    def __add__(self, other):
        result = Estimator()
        if type(other) == Estimator:
            result.value = self.value + other.value
            result.error = np.sqrt(self.error**2 + other.error**2)
        else:
            # Assume other is an int or float.
            result.value = self.value + other
            result.error = self.error
        return result
    def __radd__(self, other):
        return self.__add__(other)
    def __sub__(self, other):
        result = Estimator()
        if type(other) == Estimator:
            result.value = self.value - other.value
            result.error = np.sqrt(self.error**2 + other.error**2)
        else:
            # Assume other is an int or float.
            result.value = self.value - other
            result.error = self.error
        return result
    def __rsub__(self, other):
        result = self.__sub__(other)
        if type(other) != Estimator:
            result.value *= -1
        return result
    def __mul__(self, other):
        result = Estimator()
        if type(other) == Estimator:
            result.value = self.value*other.value
            result.error = (np.abs(result.value)
                            *np.sqrt((self.error/self.value)**2 +
                                     (other.error/other.value)**2))
        else:
            result.value = self.value*other
            result.error = self.error*other
        return result
    def __rmul__(self, other):
        return self.__mul__(other)
    def __div__(self, other):
        return self.__truediv__(other)
    def __truediv__(self, other):
        result = Estimator()
        if type(other) == Estimator:
            result.value = self.value/other.value
            result.error = (np.abs(result.value)
                            *np.sqrt((self.error/self.value)**2 +
                                     (other.error/other.value)**2))
        else:
            result.value = self.value/other
            result.error = self.error/other
        return result
    def set_format_str(self, format_str):
        self._format_str = format_str
    def __repr__(self):
        return self.__str__()
    def __str__(self):
        if self._format_str == None:
            return "%s +/- %s" % (self.value, self.error)
        return ' +/- '.join((self._format_str.format(self.value),
                             self._format_str.format(self.error)))

