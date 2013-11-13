"""
  Routine derived from James Chiang's pair_stats.py;
  Adapted to Augustin Guyonnet's algorithm of 
  pixel-correlation estimate from the diff. of
  a flat-fields pair.
  Author : Sylvain Baumont

  Use : $ python pair_pixelcovar.py 
     => test on simulated images
        $ python pair_pixelcovar.py flat1.fits flat2.fits
     => test on existing image files
"""

import copy
import numpy as np
import numpy.random as random
import pyfits

import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath

from image_utils import bias, overscan

from math import sqrt # from afw.Math ?
import os,sys

keepmask = False
#keepmask = True # Uncomment this if you wish to examine the mask images (one per ampli.)

varclip = lambda x : afwMath.makeStatistics(x, afwMath.VARIANCECLIP).getValue()

def exptime(infile):
    try:
        md = afwImage.readMetadata(infile, 1)
        return md.get('EXPTIME')
    except:
        foo = pyfits.open(infile)
        return foo[0].header['EXPTIME']


class PairCovar(object):
    def __init__(self, bias_mean, bias_stddev, flat_mean, flat_var, 
                 pixwid, covMat):
        self.bias_mean = bias_mean
        self.bias_stddev = bias_stddev
        self.flat_mean = flat_mean
        self.flat_var = flat_var
        self.pixwid = pixwid
        self.covMat = covMat
    def header(self):
        return " Ampli   Bias Mean   Bias RMS    Flat Mean   Flat Var   Covar00     Correl1X    Correl1Y    Correl%dX    Correl%dY\n ------ ----------- ----------- ----------- ----------- ----------- ----------- ----------- ----------- -----------"%(self.pixwid,self.pixwid)
    def summary(self, i=0):
        return ("%5d" % i)  + self.__repr__()
    def __repr__(self):
        format = 5*"%12.3F"+4*"%12.5F%%"
        indx = self.pixwid
        return format % (self.bias_mean, self.bias_stddev, 
                         self.flat_mean, self.flat_var,
                         self.covMat[0,0], 
                    100.*self.covMat[0,1]/(self.covMat[0,0]), 
                    100.*self.covMat[1,0]/(self.covMat[0,0]),
                    100.*self.covMat[0,indx]/(self.covMat[0,0]), 
                    100.*self.covMat[indx,0]/(self.covMat[0,0]))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def pair_covar(file1, file2, hdu=2, Nsigma=5, 
               pixwid=3, boxwid=250, boxhgt=500):
    if exptime(file1) != exptime(file2):
        raise RuntimeError("Exposure times for files %s, %s do not match"
                           % (file1, file2))

    im1 = afwImage.ImageF(file1, hdu)
    im2 = afwImage.ImageF(file2, hdu)

    if( boxwid > 0.9 * im1.getWidth()):
        boxwid = int( 0.9 * im1.getWidth() ) # avoid img. borders
    if( boxhgt > 0.95 * im1.getHeight()):
        boxhgt = int( 0.95 * im1.getHeight() )

    x0 = (im1.getWidth() - boxwid)/2 # analysis box centered on img.
    y0 = (im1.getHeight()- boxhgt)/2 

    flat_region = afwGeom.Box2I(afwGeom.Point2I(x0, y0), 
                                afwGeom.Extent2I(boxwid, boxhgt))
    # Use overscan region for bias regions.  Make a deep copy since these
    # will be used after im1 and im2 are manipulated.
    b1 = copy.deepcopy(im1.Factory(im1, overscan).getArray())
    b2 = copy.deepcopy(im2.Factory(im2, overscan).getArray())
    bmean1 = np.mean(b1)
    bmean2 = np.mean(b2)
    f1 = im1.Factory(im1, flat_region).getArray() - bmean1
    f2 = im2.Factory(im2, flat_region).getArray() - bmean2
# fratio = np.mean(f1/f2); f2 *= fratio # Better not to rescale.
    fmean = (np.mean(f1) + np.mean(f2))/2.
    fdiff = f1 - f2
    diffmed  = np.median(fdiff)
    #
    # Use clipped variance to handle bright pixels and columns
    diffvar = varclip(np.array(fdiff.flat, dtype=np.float))
    bias_rms = np.std(b1)
    #
    # Allocate the Covariance matrix :
    covMat = np.zeros((pixwid+1, pixwid+1), dtype=np.float)
    #
    # Prepare the mask of deviant pixels :
    mask = np.ones((boxhgt, boxwid), dtype=np.int)
    thresh = Nsigma*sqrt(diffvar)
    for j in range(boxhgt):
        for i in range(boxwid):
            if abs( fdiff[j,i]-diffmed ) > thresh :
                mask[j,i] = 0
    #
    # Do the job : compute the pixel-covariance as fct. of pixel offset (k,l)
    for k in range(pixwid+1):
        for l in range(pixwid+1):
            sum1 = 0.; sum2 = 0.; sum12 = 0.
            npixused = 0
            for j in range(boxhgt-l):
                for i in range(boxwid-k):
                    if( mask[j,i]==1 and mask[j+l,i+k]==1 ):
                        val1 = fdiff[j,i]; val2 = fdiff[j+l,i+k]
                        sum1 += val1; sum2 += val2
                        sum12 += val1*val2
                        npixused += 1
            
            covMat[l,k] = (sum12 - sum1*sum2/npixused)/npixused

    if keepmask :
        outfile = "mask_%dsig_hdu%d.fits"%(Nsigma, hdu-1)
        output = pyfits.HDUList()
        output.append(pyfits.PrimaryHDU(mask))
        output.writeto(outfile, clobber=True)

    return PairCovar(bmean1, bias_rms, fmean, diffvar/2., pixwid, covMat)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if __name__ == '__main__':

    if( len(sys.argv)>2 ):
        file1 = sys.argv[1]
        file2 = sys.argv[2]
        nHdus = 16
        print " Use files '%s' & '%s'"%(file1,file2)
    else:
        from simulation.sim_tools import simulateFlat

        file1 = 'test_flat1.fits'
        file2 = 'test_flat2.fits'

        nHdus = 4
        simulateFlat(file1, 15000, 5, hdus=nHdus)
        simulateFlat(file2, 15000, 5, hdus=nHdus)

    for hdu in range(nHdus):
        my_pair_cov = pair_covar(file1, file2, hdu=hdu+2)
        if hdu == 0:
            print my_pair_cov.header()
        print my_pair_cov.summary( hdu+1 )
