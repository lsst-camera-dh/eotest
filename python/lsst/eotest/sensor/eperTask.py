#/usr/bin/env python

# Charge transfer efficiency by EPER, now as a pipe task!

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

import sys
import numpy as np
import argparse
from MaskedCCD import MaskedCCD
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath

import lsst.eotest.image_utils as imutils
from lsst.eotest.Estimator import Estimator
from AmplifierGeometry import makeAmplifierGeometry

class SubImage(object):
    """Functor to produce sub-images depending on scan direction."""
    def __init__(self, ccd, amp, overscans, task):
        geom = ccd.amp_geom
        self.ccd = ccd
        self.imaging = geom.imaging
        self.image = ccd[amp]   # This is the masked image for the desired amp.
        if task.config.direction == 'p':
            self._bbox = self._parallel_box
            llc = afwGeom.Point2I(geom.parallel_overscan.getMinX(),
                                  geom.parallel_overscan.getMinY() + overscans)
            urc = geom.parallel_overscan.getCorners()[2]
            self._bias_reg = afwGeom.Box2I(llc, urc)
            self.lastpix = self.imaging.getMaxY()
        elif task.config.direction == 's':
            self._bbox = self._serial_box
            llc = afwGeom.Point2I(geom.serial_overscan.getMinX() + overscans,
                                  geom.serial_overscan.getMinY())
            urc = geom.serial_overscan.getCorners()[2]
            #
            # Omit the last 4 columns to avoid the bright column in the
            # last overscan column in the e2v vendor data.
            #
            urc[0] -= 4
            self._bias_reg = afwGeom.Box2I(llc, urc)
            self.lastpix = self.imaging.getMaxX()
        else:
            task.log.error("Unknown scan direction: " + str(direction))
            sys.exit(1)
    def bias_est(self, statistic=afwMath.MEAN, gain=1):
        subim = self.image.Factory(self.image, self._bias_reg)
        # Set bias region error to zero since it isn't governed by
        # Poisson statistics and cannot be estimated only from a
        # medianed, bias-subtracted superflat frame stack.
        bias_estimate = Estimator()
        bias_estimate.value = \
            gain*afwMath.makeStatistics(subim, statistic).getValue()
        bias_estimate.error = \
            gain*afwMath.makeStatistics(subim, afwMath.STDEV).getValue()
        return bias_estimate
    def __call__(self, start, end=None):
        if end is None:
            end = start
        my_exp = self.image.Factory(self.image, self._bbox(start, end))
        return my_exp
    def _parallel_box(self, start, end):
        llc = afwGeom.PointI(self.imaging.getMinX(), start)
        urc = afwGeom.PointI(self.imaging.getMaxX(), end)
        return afwGeom.BoxI(llc, urc)
    def _serial_box(self, start, end):
        llc = afwGeom.PointI(start, self.imaging.getMinY())
        urc = afwGeom.PointI(end, self.imaging.getMaxY())
        return afwGeom.BoxI(llc, urc)

class EPERConfig(pexConfig.Config):
    """Configuration for the EPERTask."""
    direction = pexConfig.Field("Select either parallel or serial direction",
                                str, default="p")
    verbose = pexConfig.Field("Turn verbosity on", bool, default=True)
    cti = pexConfig.Field('Return CTI instead of CTE', bool, default=False)

class EPERTask(pipeBase.Task):
    """Task to calculate either parallel or serial charge transfer
       efficiency via EPER."""
    ConfigClass = EPERConfig
    _DefaultName = "eper"

    @pipeBase.timeMethod
    def run(self, infilename, amps, overscans, gains=None):
        if not infilename:
            self.log.error("Please specify an input file path.")
            sys.exit(1)
        if gains is None:
            gains = dict([(amp, 1) for amp in amps])

        ccd = MaskedCCD(infilename)
        # iterate through amps
        cte = {}
        bias_estimates = {}
        for amp in amps:
            subimage = SubImage(ccd, amp, overscans, self)
            lastpix = subimage.lastpix

            # find signal in last image vector (i.e., row or column)
            last_im = Estimator(subimage(lastpix), ccd.stat_ctrl,
                                gain=gains[amp])
            if self.config.verbose:
                self.log.info("Last imaging row/column = " + str(last_im))

            # find signal in each overscan vector
            overscan_ests = []
            for i in range(1, overscans+1):
                overscan_ests.append(Estimator(subimage(lastpix+i),
                                               ccd.stat_ctrl, gain=gains[amp]))
            if self.config.verbose:
                self.log.info("Overscan values = " + str(overscan_ests))

            # sum medians of first n overscan rows
            summed = sum(overscan_ests)
            if self.config.verbose:
                self.log.info("summed overscans = " + str(summed))

            # Find bias level, use afwMath.MEANCLIP to avoid
            # contribution of bad pixels in overscan regions (e.g.,
            # e2v vendor data).
            bias_est = subimage.bias_est(gain=gains[amp],
                                         statistic=afwMath.MEANCLIP)
            bias_estimates[amp] = bias_est
            if self.config.verbose:
                self.log.info("bias value = " + str(bias_est))

            # signal = last - bias
            sig = last_im - bias_est

            # trailed = sum(last2) - bias
            trailed = summed - overscans*bias_est

            # charge loss per transfer = (trailed/signal)/N
            chargelosspt = (trailed/sig)/(lastpix + 1.)

            if self.config.cti:
                cte[amp] = chargelosspt
                cte[amp].set_format_str("{0:.5e}")
            else:
                cte[amp] = 1. - chargelosspt
                cte[amp].set_format_str("{0:.16f}")
            if self.config.verbose:
                if self.config.cti:
                    self.log.info('cti, amp ' + str(amp) + " = "
                                  + str(cte[amp]) + '\n')
                else:
                    self.log.info('cte, amp ' + str(amp) + " = "
                                  + str(cte[amp]) + '\n')
        return cte, bias_estimates

if __name__ == '__main__':
    #import pdb; pdb.set_trace()
    parser = argparse.ArgumentParser(description='Calculate either parallel or serial CTE via EPER.')
    parser.add_argument('infilename', help="image file to be used for analysis")
    parser.add_argument('-o', '--overscans', help = "number of overscan rows/columns to use", type=int, default=3)
    parser.add_argument('-d', '--direction', help="specify either parallel ('p') or serial ('s') direction", default='p')
    parser.add_argument('-a', '--amps', help="amps to be analyzed, separated by a space", type=int, nargs = '+', default=range(1, 17))
    parser.add_argument('-v', '--verbose', help="turn verbosity on", action='store_true', default=False)
    parser.add_argument('-i', '--cti', help='return CTI (not CTE)',
                        action='store_true', default=False)
    args = parser.parse_args()

    task = EPERTask()
    task.config.direction = args.direction
    task.config.verbose = args.verbose
    task.config.cti = args.cti

    task.run(args.infilename, args.amps, args.overscans)
