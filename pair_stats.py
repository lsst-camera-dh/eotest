import copy
import numpy as np
import numpy.random as random
import pyfits

import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage

from image_utils import bias, overscan

def exptime(infile):
    try:
        md = afwImage.readMetadata(infile, 1)
        return md.get('EXPTIME')
    except:
        foo = pyfits.open(infile)
        return foo[0].header['EXPTIME']

class PairStats(object):
    def __init__(self, bias_mean, bias_stddev, flat_mean, flat_var, 
                 gain, noise):
        self.bias_mean = bias_mean
        self.bias_stddev = bias_stddev
        self.flat_mean = flat_mean
        self.flat_var = flat_var
        self.gain = gain
        self.noise = noise
    def header(self):
        return """ Exp     Bias Mean   Bias RMS    Flat Mean   Flat Var   Gain e-/DN   Noise e-  
 ------ ----------- ----------- ----------- ----------- ----------- -----------"""
    def summary(self, i=0):
        return ("%5.2f" % i)  + self.__repr__()
    def __repr__(self):
        format = 6*"%12.3F"
        return format % (self.bias_mean, self.bias_stddev, 
                         self.flat_mean, self.flat_var,
                         self.gain, self.noise)

def pair_stats(file1, file2, hdu=2):
    if exptime(file1) != exptime(file2):
        raise RuntimeError("Exposure times for files %s, %s do not match"
                           % (file1, file2))
    flat_region = afwGeom.Box2I(afwGeom.Point2I(200, 900), 
                                afwGeom.Extent2I(100, 100))
    im1 = afwImage.ImageF(file1, hdu)
    im2 = afwImage.ImageF(file2, hdu)
    # Use overscan region for bias regions.  Make a deep copy since these
    # will be used after im1 and im2 are manipulated.
    b1 = copy.deepcopy(im1.Factory(im1, overscan).getArray())
    b2 = copy.deepcopy(im2.Factory(im2, overscan).getArray())
    bmean = (np.mean(b1) + np.mean(b2))/2.
    f1 = im1.Factory(im1, flat_region).getArray() - bmean
    f2 = im2.Factory(im2, flat_region).getArray() - bmean
    fratio = np.mean(f1/f2)
    f2 *= fratio
    fmean = (np.mean(f1) + np.mean(f2))/2.
    fdiff = f1 - f2
    bdiff = b1 - b2
    fvar = np.var(fdiff)/2.
    bvar = np.var(bdiff)/2.
    gain = fvar/fmean
#    gain = (fvar - bvar)/fmean  # seems like we should subtract bvar.
    bias_rms = np.std(b1)
    noise = gain*np.std(b1)
    return PairStats(bmean, bias_rms, fmean, fvar, gain, noise), b1, b2

if __name__ == '__main__':
    from simulation.sim_tools import simulateFlat

    file1 = 'test_flat1.fits'
    file2 = 'test_flat2.fits'

    simulateFlat(file1, 200, 5, hdus=16)
    simulateFlat(file2, 200, 5, hdus=16)

    for hdu in range(16):
        my_pair_stats, b1, b2 = pair_stats(file1, file2, hdu=hdu+2)
        if hdu == 0:
            print my_pair_stats.header()
        print my_pair_stats.summary()
