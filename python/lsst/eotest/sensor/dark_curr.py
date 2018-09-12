"""
@brief Re-implementation of P. Doherty's IDL function to compute
dark current.
"""
from __future__ import print_function
from builtins import zip
from builtins import range
import numpy as np
import numpy.random as random

import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom

from image_utils import fits_median, unbias_and_trim


def dark_curr(files, hdu=2, gain=1, count=1000, dx=100, dy=100, seed=None):
    random.seed(seed)
    exptime = afwImage.readMetadata(files[0], 1).get('EXPTIME')
    im = unbias_and_trim(fits_median(files, hdu=hdu))

    # Generate dx by dy boxes at random locations to perform
    # estimates, then take the median.  This avoids bright defects.
    xarr = random.randint(im.getWidth() - dx - 1, size=count)
    yarr = random.randint(im.getHeight() - dy - 1, size=count)

    signal = []
    for x, y in zip(xarr, yarr):
        bbox = afwGeom.Box2I(afwGeom.Point2I(int(x), int(y)),
                             afwGeom.Extent2I(dx, dy))
        subim = im.Factory(im, bbox)
        signal.append(np.mean(subim.getArray()))
    dark_current = np.median(signal)*gain/exptime

    return dark_current


if __name__ == '__main__':
    from simulation.sim_tools import simulateDark

    hdus = 2
    dark_current = 3
    exptime = 12
    files = []
    for i in range(4):
        outfile = 'test_dark%02i.fits' % i
        simulateDark(outfile, dark_current, exptime, hdus=hdus)
        files.append(outfile)

    for hdu in range(hdus):
        print(hdu, dark_curr(files, hdu=hdu+2, count=100))
