"""
@brief Calculate the gain and noise of a CCD camera system by
examining two flat field images. The calculation is the standard
mean/variance thing.
"""
import numpy.random as random

import lsst.afw.image as afwImage

from image_utils import bias

def flat_gain(file1, file2, hdu=2, count=1000, seed=19031, dx=100, dy=100):
    im1 = afwImage.ImageF(file1, hdu=hdu)
    im2 = afwImage.ImageF(file2, hdu=hdu)

    # Unbias using the mean bias of the two images.
    bmean = (bias(im1) + bias(im2))/2.
    im1 -= bmean
    im2 -= bmean

    # Sample detector at size=count locations.
    xarr = random.randint(im1.getWidth() - dx - 1, size=count)
    yarr = random.randint(im1.getHeight() - dy - 1, size=count)

    gains = []
    for x, y in zip(xarr, yarr):
        bbox = afwGeom.Box2I(afwGeom.Point2I(x, y), afwGeom.Extent2I(dx, dy))
        imarr1 = im1.Factory(im1, bbox).getArray()
        imarr2 = im2.Factory(im2, bbox).getArray()

        # Calculate flat ratio and correct subarray of image 2.
        fratio = np.mean(imarr1/imarr2)
        imarr2 *= fratio

        # Calculate the mean value of the flat field images.
        fmean = (np.mean(imarr1) + np.mean(imarr2))/2.
        
        # Calculate the variance of the flat difference image.
        fvar = np.std(imarr1 - imarr2)**2/2.

        gains.append(fmean/fvar)
    gain = np.median(gains)
    return gain

if __name__ == '__main__':
    import glob
    files = glob.glob('data/flat*.fits')[:2]

    seed = 1304291
    for hdu in range(2, 18):
        print hdu, flat_gain(files[0], files[1], count=1000, seed=seed)
