"""
@brief Re-implementation of P. Doherty's IDL function to compute
dark current.
"""
import numpy as np
import numpy.random as random

import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom

from image_utils import fits_median, unbias_and_trim

def dark_curr(files, hdu=2, gain=1, count=1000, 
              nx=400, ny=900, x0=10, y0=0, dx=100, dy=100, seed=101):
    random.seed(seed)
    im = unbias_and_trim(fits_median(files, hdu=hdu))

    # Generate random locations to perform estimates, then take the
    # median.  This avoids bright defects.
    xarr = random.randint(nx, size=count) + x0
    yarr = random.randint(ny, size=count) + y0

    signal = []
    for x, y in zip(xarr, yarr):
        llc = afwGeom.Point2I(int(x), int(y))
        extent = afwGeom.Extent2I(dx, dy)
        bbox = afwGeom.Box2I(llc, extent)
        subim = im.Factory(im, bbox)
        signal.append(np.mean(subim.getArray()))
    dark_current = np.median(signal)*gain

    return dark_current

if __name__ == '__main__':
    import glob
    files = glob.glob('data/dark*.fits')
    seed = 10184791
    for hdu in range(2, 6):
        print hdu, dark_curr(files, hdu=hdu, count=100, seed=seed)
