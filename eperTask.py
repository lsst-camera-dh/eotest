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

class SubImage(object):
    """Functor to produce sub-image depending on scan direction."""
    def __init__(self, imfile, amp, startimagepix, endimagepix, task):
        self.exp = afwImage.ExposureF(imfile, int(amp) + 1)
        self.startimagepix = startimagepix
        self.endimagepix = endimagepix
        self.direction = task.config.direction
        if direction not in ('p', 's'):
            task.log.error("Unknown scan direction: " + str(direction))
            sys.exit(1)
    def __call__(self, start, end=None):
        if end is None:
            end = start + 1
        my_exp = self.exp.Factory(self.exp, self._bbox(start, end))
        return my_exp.getMaskedImage().getImage()
    def _bbox(self, start, end):
        if self.direction == 'p':
            llc = afwGeom.PointI(self.startimagepix, start)
            urc = afwGeom.PointI(self.endimagepix, end)
        else:
            llc = afwGeom.PointI(start, self.startimagepix)
            urc = afwGeom.PointI(end, self.endimagepix)
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
    def run(self, infilename, amps, startimagepix, endimagepix, 
            lastpix, biasstart, biasend, overscans):
        if not infilename:
            self.log.error("Please specify an input file path.")
            sys.exit(1)

        median = lambda x : afwMath.makeStatistics(x, afwMath.MEDIAN).getValue()
		
        cte = []
        # iterate through amps
        for amp in amps:
            subimage = SubImage(infilename, amp, startimagepix,
                                endimagepix, self)
            # find signal in last image row/column
            last = subimage(lastpix)
            lastmed = median(last)
			
            if self.config.verbose:
                print "lastmed = " + str(lastmed)
		
            # find median signal in each overscan row
            overscanmeds = np.zeros(overscans)
            for i in range(1, overscans+1):
                eachover = subimage(lastpix + i)
                overscanmeds[i-1] = median(eachover)
			
            if self.config.verbose:
                print "Overscan medians = " + str(overscanmeds)
		
            # sum medians of first n overscan rows
            summed = np.sum(overscanmeds)
			
            if self.config.verbose:
                print "summed = " + str(summed)

            # find signal in bias
            bias = subimage(biasstart, biasend)
            biasmed = median(bias)
			
            if self.config.verbose:
                print "biasmed = " + str(biasmed)

            # signal = last -bias
            sig = lastmed - biasmed

            # trailed = sum(last2) - bias
            trailed = summed - 1.0*overscans*biasmed		

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
    parser.add_argument('-l', '--last', help="last row/column of image", type=int, default = 2000)
    parser.add_argument('-s', '--startimagepix', help = "first column/row of image", type=int, default = 10)
    parser.add_argument('-e', '--endimagepix', help = "last column/row of image", type=int, default = 522)
    parser.add_argument('-b', '--biasstart', help = "first row/column of bias region", type=int, default = 2010)
    parser.add_argument('-i', '--biasend',  help = "last row/column of bias region", type=int, default = 2015)
    parser.add_argument('-o', '--overscans', help = "number of overscan rows/columns to use", type=int, default = 3)
    parser.add_argument('-d', '--direction', help="specify either parallel ('p') or serial ('s') direction", default = 'p')
    parser.add_argument('-a', '--amps', help="amps to be analyzed, separated by a space", type=int, nargs = '+', default = range(1,17))
    parser.add_argument('-v', '--verbose', help="turn verbosity on", action='store_true', default=False)
    args = parser.parse_args()
	
    task = EPERTask()
    task.config.direction = args.direction
    task.config.verbose = args.verbose
    task.run(args.infilename, args.amps, args.startimagepix, 
             args.endimagepix, args.last, args.biasstart, args.biasend, 
             args.overscans)
