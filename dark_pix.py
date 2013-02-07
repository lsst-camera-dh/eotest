import lsst.pipe.base as pipeBase

import sys
import numpy as np
import argparse
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath
import lsst.afw.detection as afwDetect
import lsst.afw.display.ds9 as ds9


def dark_pix(infilename, percent, amps):
    """ List pixels with counts less than a specified percentage of the
        median flux, for the specified amps. """

    for amp in amps:
        #read in and trim image area
        im = afwImage.ExposureF(infilename, amp+1, \
                      afwGeom.BoxI(afwGeom.PointI(10, 0), \
                      afwGeom.PointI(521, 2001))).getMaskedImage().getImage()
        
        #find median of image
        median = afwMath.makeStatistics(im, afwMath.MEDIAN).getValue()
        thresh = median*percent/100.0
        #print thresh
        
        #find pixels less than _ percent of median
        imarr = im.getArray()
        darkpix = np.where(imarr <= thresh)
        
#        p.figure()
#        p.plot(darkpix)
        
        #turn x,y into a list of pixel coords
        pixlist = np.transpose(darkpix)

        
        print thresh
        print pixlist
        print len(pixlist)
        
        return pixlist

if __name__ == '__main__':
    parser = argparse.ArgumentParser(\
                      description='Find the locations of dark pixels.')
    parser.add_argument('-i', '--infile', help="path to input image file" \
                        type=str)
    parser.add_argument('-o', '--outfile', help="output results text file" \
                        type=str)
    parser.add_argument('-p', '--percent', \
                        help="Percentage of median to use as threshold", \
                        type=float, default = 80.0)
    parser.add_argument('-a', '--amps', \
                        help="amps to be analyzed, separated by a space", \
                        type=int, nargs = '+', default = range(1,17))
    parser.add_argument('-v', '--verbose', help="turn verbosity on",\
                        action='store_true', default=False)
    args = parser.parse_args()

    
    if args.outfile:
        outfile = open(args.outfile, "w")
    else:
        outfile = None

    im = dark_pix(args.infilename, args.percent, args.amps)

    if outfile:
        outfile.close()

