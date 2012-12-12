"""
@brief Re-implementation of P. Doherty's IDL function to compute
dark current.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import numpy as num
import numpy.random as random

import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom

from fits_median import fits_median
from sensorImage import SensorImage

def dark_curr(files, hdu=2, gain=1, count=1000, 
              nx=400, ny=900, x0=10, y0=0, dx=100, dy=100, seed=101):
    random.seed(seed)
    im = SensorImage(fits_median(files, hdu=hdu)).unbias()

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
        signal.append(num.mean(subim.getArray()))
    dark_current = num.median(signal)*gain

    return dark_current

if __name__ == '__main__':
    import glob
    files = glob.glob('data/dark*.fits')
    seed = 10184791
    for hdu in range(2, 6):
        print hdu, dark_curr(files, hdu=hdu, count=100, seed=seed)
