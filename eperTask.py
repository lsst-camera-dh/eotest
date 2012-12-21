#/usr/bin/env python

# Charge transfer efficiency by EPER, now as a pipe task!

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

import sys
import numpy as np
import argparse
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath

from image_utils import parallel_overscan, serial_overscan, imaging

median = lambda x : afwMath.makeStatistics(x, afwMath.MEDIAN).getValue()

class SubImage(object):
    """Functor to produce sub-images depending on scan direction."""
    def __init__(self, imfile, amp, overscans, task):
        self.exp = afwImage.ExposureF(imfile, int(amp) + 1)
        if task.config.direction == 'p':
            self._bbox = self._parallel_box
            llc = afwGeom.Point2I(parallel_overscan.getMinX(),
                                  parallel_overscan.getMinY() + overscans)
            urc = parallel_overscan.getCorners()[2]
            self._bias_reg = afwGeom.Box2I(llc, urc)
            self.lastpix = imaging.getMaxY()
        elif task.config.direction == 's':
            self._bbox = self._serial_box
            llc = afwGeom.Point2I(serial_overscan.getMinX() + overscans,
                                  serial_overscan.getMinY())
            urc = serial_overscan.getCorners()[2]
            self._bias_reg = afwGeom.Box2I(llc, urc)
            self.lastpix = imaging.getMaxX()
        else:
            task.log.error("Unknown scan direction: " + str(direction))
            sys.exit(1)
    def bias_med(self):
        subim = self.exp.Factory(self.exp, self._bias_reg)
        bias = subim.getMaskedImage().getImage()
        return median(bias)
    def __call__(self, start, end=None):
        if end is None:
            end = start + 1
        my_exp = self.exp.Factory(self.exp, self._bbox(start, end))
        return my_exp.getMaskedImage().getImage()
    def _parallel_box(self, start, end):
        llc = afwGeom.PointI(imaging.getMinX(), start)
        urc = afwGeom.PointI(imaging.getMaxX(), end)
        return afwGeom.BoxI(llc, urc)
    def _serial_box(self, start, end):
        llc = afwGeom.PointI(start, imaging.getMinY())
        urc = afwGeom.PointI(end, imaging.getMaxY())
        return afwGeom.BoxI(llc, urc)

class EPERConfig(pexConfig.Config):
    """Configuration for the EPERTask."""
    direction = pexConfig.Field("Select either parallel or serial direction", 
                                str, default="p")
    verbose = pexConfig.Field("Turn verbosity on", bool, default=True)

class EPERTask(pipeBase.Task):
    """Task to calculate either parallel or serial charge transfer 
       efficiency via EPER."""
    ConfigClass = EPERConfig
    _DefaultName = "eper"
	 
    @pipeBase.timeMethod
    def run(self, infilename, amps, overscans):
        if not infilename:
            self.log.error("Please specify an input file path.")
            sys.exit(1)

        # iterate through amps
        cte = []
        for amp in amps:
            subimage = SubImage(infilename, amp, overscans, self)
            lastpix = subimage.lastpix

            # find signal in last image row/column
            lastmed = median(subimage(lastpix))
            if self.config.verbose:
                print "lastmed = " + str(lastmed)
		
            # find median signal in each overscan row
            overscanmeds = np.zeros(overscans)
            for i in range(1, overscans+1):
                overscanmeds[i-1] = median(subimage(lastpix + i))
            if self.config.verbose:
                print "Overscan medians = " + str(overscanmeds)
		
            # sum medians of first n overscan rows
            summed = np.sum(overscanmeds)
            if self.config.verbose:
                print "summed = " + str(summed)

            # find signal in bias
            biasmed = subimage.bias_med()
            if self.config.verbose:
                print "biasmed = " + str(biasmed)

            # signal = last - bias
            sig = lastmed - biasmed

            # trailed = sum(last2) - bias
            trailed = summed - overscans*biasmed		

            # charge loss per transfer = (trailed/signal)/N
            chargelosspt = (trailed/sig)/(lastpix + 1.)

            cte.append(1. - chargelosspt)
            if self.config.verbose:
                print 'cte, amp ' + str(amp) + " = " + '{0:.16f}'.format(cte[-1])
        return cte
			
if __name__ == '__main__':
    #import pdb; pdb.set_trace()
    parser = argparse.ArgumentParser(description='Calculate either parallel or serial CTE via EPER.')
    parser.add_argument('infilename', help="image file to be used for analysis")
    parser.add_argument('-o', '--overscans', help = "number of overscan rows/columns to use", type=int, default=3)
    parser.add_argument('-d', '--direction', help="specify either parallel ('p') or serial ('s') direction", default='p')
    parser.add_argument('-a', '--amps', help="amps to be analyzed, separated by a space", type=int, nargs = '+', default=range(1, 17))
    parser.add_argument('-v', '--verbose', help="turn verbosity on", action='store_true', default=False)
    args = parser.parse_args()
	
    task = EPERTask()
    task.config.direction = args.direction
    task.config.verbose = args.verbose

    task.run(args.infilename, args.amps, args.overscans)
