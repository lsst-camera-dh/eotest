"""
@brief Calculate the gain and noise of a CCD camera system by
examining two flat field images. The calculation is the standard
mean/variance thing.
"""
import numpy as np
import numpy.random as random
import lsst.eotest.image_utils as imutils
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.pex.exceptions as pexExcept

def flat_gain(image1, image2, count=1000, dx=100, dy=100, binsize=1,
              seed=None):
    """
    Calculate the gain and noise of a CCD camera system by
    examining two flat field images. The calculation is the standard
    mean/variance thing.
    """
    # If seed is None, the seed is generated from /dev/urandom.
    random.seed(seed)

    # Unbias using the mean bias of the two images.
    bmean = (imutils.bias(image1) + imutils.bias(image2))/2.
    image1 -= bmean
    image2 -= bmean

    # Trim prescan and overscan.
    image1 = imutils.trim(image1)
    image2 = imutils.trim(image2)

    # Rebin into binsize x binsize pixels.
    im1 = imutils.rebin(image1, binsize)
    im2 = imutils.rebin(image2, binsize)

    if dx > im1.getWidth():
        dx = im1.getWidth()
    if dy > im1.getHeight():
        dy = im1.getHeight()

    # Sample detector at size=count locations.
    try:
        xarr = random.randint(im1.getWidth() - dx - 1, size=count)
    except ValueError:
        # Rebinned image width smaller than requested dx, so just
        # generate subregions using the full x extent.
        xarr = np.zeros(count, dtype=np.int)
    try:
        yarr = random.randint(im1.getHeight() - dy - 1, size=count)
    except ValueError:
        # Rebinned image height smaller than requested dy, so just
        # generate subregions using the full y extent.
        yarr = np.zeros(count, dtype=np.int)

    gains = []
    ntrial = 0
    exception_count = 0
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
    gain = 1./np.median(gains)  # gain in Ne/DN
    return gain, im1, im2

if __name__ == '__main__':
    import os
    from sim_tools import simulateFlat

    file1 = 'test_flat1.fits'
    file2 = 'test_flat2.fits'

    hdus = 2
    simulateFlat(file1, 4000, 5, hdus=hdus)
    simulateFlat(file2, 4000, 5, hdus=hdus)

    for amp in range(1, hdus+1):
        image1 = afwImage.ImageF(file1, imutils.dm_hdu(amp))
        image2 = afwImage.ImageF(file2, imutils.dm_hdu(amp))
        gain, im1, im2 = flat_gain(image1, image2, count=1000)
        print amp, gain

    os.remove(file1)
    os.remove(file2)
