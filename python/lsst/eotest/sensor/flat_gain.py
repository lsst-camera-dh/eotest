"""
@brief Calculate the gain and noise of a CCD camera system by
examining two flat field images. The calculation is the standard
mean/variance thing.
"""
import numpy as np
import numpy.random as random

import lsst.eotest.image_utils as imutils
import lsst.eotest.utilLib as testUtils

import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage

def flat_gain(file1, file2, hdu=2, count=1000, seed=None, dx=100, dy=100,
              binsize=3):
    """
    Calculate the gain and noise of a CCD camera system by
    examining two flat field images. The calculation is the standard
    mean/variance thing.
    """
    # If seed is None, the seed is generated from /dev/urandom.
    random.seed(seed)

    dx /= binsize
    dy /= binsize

    im1 = afwImage.ImageF(file1, hdu)
    im2 = afwImage.ImageF(file2, hdu)

    # Unbias using the mean bias of the two images.
    bmean = (imutils.bias(im1) + imutils.bias(im2))/2.
    im1 -= bmean
    im2 -= bmean

    # Trim prescan and overscan.
    im1 = imutils.trim(im1)
    im2 = imutils.trim(im2)

    # Rebin into binsize x binsize pixels.
    im1 = testUtils.ImageTools_rebin(im1, binsize)
    im2 = testUtils.ImageTools_rebin(im2, binsize)

    # Sample detector at size=count locations.
    xarr = random.randint(im1.getWidth() - dx - 1, size=count)
    yarr = random.randint(im1.getHeight() - dy - 1, size=count)

    gains = []
    ntrial = 0
    for x, y in zip(xarr, yarr):
        bbox = afwGeom.Box2I(afwGeom.Point2I(int(x), int(y)), 
                             afwGeom.Extent2I(dx, dy))
        imarr1 = im1.Factory(im1, bbox).getArray()
        imarr2 = im2.Factory(im2, bbox).getArray()
        # Calculate flat ratio and correct subarray of image 2.
        fratio = np.mean(imarr1/imarr2)
        imarr2 *= fratio
        # Calculate the mean value of the flat field images.
        fmean = (np.mean(imarr1) + np.mean(imarr2))/2.
        # Calculate the variance of the flat difference image.
        fvar = np.var(imarr1 - imarr2)/2.
        gains.append(fvar/fmean)
    gain = np.median(gains)
    return gain, im1, im2

if __name__ == '__main__':
    from simulation.sim_tools import simulateFlat
    
    file1 = 'test_flat1.fits'
    file2 = 'test_flat2.fits'

    hdus = 2
    simulateFlat(file1, 40, 5, hdus=hdus)
    simulateFlat(file2, 40, 5, hdus=hdus)

    for hdu in range(hdus):
        print hdu, flat_gain(file1, file2, hdu=hdu+2, count=1000)
