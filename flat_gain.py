"""
@brief Calculate the gain and noise of a CCD camera system by
examining two flat field images. The calculation is the standard
mean/variance thing.
"""
import numpy as np
import numpy.random as random

import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage

from image_utils import bias, trim

def flat_gain(file1, file2, hdu=2, count=1000, seed=None, dx=100, dy=100):
    """
    Calculate the gain and noise of a CCD camera system by
    examining two flat field images. The calculation is the standard
    mean/variance thing.
    """
    # If seed is None, the seed is generated from /dev/urandom.
    random.seed(seed)

    im1 = afwImage.ImageF(file1, hdu)
    im2 = afwImage.ImageF(file2, hdu)

    # Unbias using the mean bias of the two images.
    bmean = (bias(im1) + bias(im2))/2.
    im1 -= bmean
    im2 -= bmean

    # Trim prescan and overscan.
    im1 = trim(im1)
    im2 = trim(im2)

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
    return gain

if __name__ == '__main__':
    from simulation.sim_tools import simulateFlat
    
    file1 = 'test_flat1.fits'
    file2 = 'test_flat2.fits'

    hdus = 2
    simulateFlat(file1, 40, 5, hdus=hdus)
    simulateFlat(file2, 40, 5, hdus=hdus)

    for hdu in range(hdus):
        print hdu, flat_gain(file1, file2, hdu=hdu+2, count=1000)
