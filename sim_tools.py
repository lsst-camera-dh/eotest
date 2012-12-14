import numpy as np
import numpy.random as random

import pyfits

import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.display.ds9 as ds9

from image_utils import imaging

ccd_segment = afwGeom.Box2I(afwGeom.Point2I(0, 0), afwGeom.Extent2I(543, 2023))

class SegmentExposure(object):
    def __init__(self, exptime=1, bbox=ccd_segment):
        self.exptime = exptime
        self.image = afwImage.ImageF(bbox)
        self.imarr = self.image.Factory(self.image, imaging).getArray()
        self.ny, self.nx = self.imarr.shape
        self.npix = self.nx*self.ny
    def add_bias(self, level=1e4, sigma=10):
        fullarr = self.image.getArray()
        ny, nx = fullarr.shape
        bias_arr = np.array(random.normal(level, sigma, nx*ny),
                            dtype=np.int).reshape(ny, nx)
        fullarr += bias_arr
    def add_dark_current(self, level):
        dark_arr = self._poisson_imarr(level*self.exptime)
        self.imarr += dark_arr
    def expose_flat(self, level, gain):
        flat_arr = self._poisson_imarr(level*self.exptime)*gain
        self.imarr += flat_arr
    def _poisson_imarr(self, level):
        return random.poisson(level, self.npix).reshape(self.ny, self.nx)

def writeFits(ccd_segments, outfile, clobber=True):
    output = pyfits.HDUList()
    output.append(pyfits.PrimaryHDU())
    output[0].header["EXPTIME"] = ccd_segments[0].exptime
    for segment in ccd_segments:
        output.append(pyfits.ImageHDU(data=segment.image.getArray()))
    output.writeto(outfile, clobber=clobber)
    
def simulateDark(outfile, dark_curr, exptime=1, hdus=16, verbose=True):
    if verbose:
        print "simulating dark:", outfile
    segments = []
    for i in range(hdus):
        if verbose:
            print "HDU", i
        seg = SegmentExposure(exptime)
        seg.add_bias()
        seg.add_dark_current(dark_curr)
        segments.append(seg)
    writeFits(segments, outfile)

def simulateFlat(outfile, level, gain, dark_curr=1, exptime=1, hdus=16,
                 verbose=True):
    if verbose:
        print "simulating flat:", outfile
    segments = []
    for i in range(hdus):
        if verbose:
            print "HDU", i
        seg = SegmentExposure(exptime)
#        seg.add_bias()
#        seg.add_dark_current(dark_curr)
        seg.expose_flat(level, gain)
        segments.append(seg)
    writeFits(segments, outfile)

if __name__ == '__main__':
    seg = SegmentExposure()
    
    seg.add_bias(1e4, 10)
    seg.add_dark_current(300)
    seg.expose_flat(200, 5)

    writeFits((seg,), 'test_image.fits')

#    ds9.mtv(exp.image)
