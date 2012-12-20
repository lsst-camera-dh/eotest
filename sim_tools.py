"""
@brief Tools to create simulated CCD segment exposures under ideal
conditions.  Darks, flats, Fe55, etc..
"""
import os

import numpy as np
import numpy.random as random

import pyfits

import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.display.ds9 as ds9

from image_utils import imaging, full_segment

class SegmentExposure(object):
    def __init__(self, exptime=1, bbox=full_segment):
        self.exptime = exptime
        self.image = afwImage.ImageF(bbox)
        self.imarr = self.image.Factory(self.image, imaging).getArray()
        self.ny, self.nx = self.imarr.shape
        self.npix = self.nx*self.ny
        self._sigma = -1
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
    def sigma(self):
        if self._sigma == -1:
            self._sigma = np.std(self.imarr)
        return self._sigma
    def add_bright_cols(self, ncols=1, nsig=5):
        bright_cols = np.arange(self.nx)
        random.shuffle(bright_cols)
        bright_cols = bright_cols[:ncols]
        for i in bright_cols:
            self.imarr[:, i] += nsig*self.sigma()
        return bright_cols
    def add_bright_pix(self, npix=100, nsig=5):
        bright_pix = np.concatenate((np.ones(npix), np.zeros(self.npix-npix)))
        random.shuffle(bright_pix)
        bright_pix = bright_pix.reshape(self.ny, self.nx)
        self.imarr += bright_pix*nsig*self.sigma()
        return np.where(bright_pix == 1)

def writeFits(ccd_segments, outfile, clobber=True):
    output = pyfits.HDUList()
    output.append(pyfits.PrimaryHDU())
    output[0].header["EXPTIME"] = ccd_segments[0].exptime
    for segment in ccd_segments:
        output.append(pyfits.ImageHDU(data=segment.image.getArray()))
    if clobber:
        try:
            os.remove(outfile)
        except OSError:
            pass
    output.writeto(outfile)
    
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
        seg.add_bias()
        seg.add_dark_current(dark_curr)
        seg.expose_flat(level, gain)
        segments.append(seg)
    writeFits(segments, outfile)

if __name__ == '__main__':
    seg = SegmentExposure()
    
    seg.add_bias(1e4, 10)
    seg.add_dark_current(300)
    seg.expose_flat(200, 5)
    cols = seg.add_bright_cols(ncols=1, nsig=5)
    pix = seg.add_bright_pix(npix=100, nsig=10)

    writeFits((seg,), 'test_image.fits')

#    ds9.mtv(exp.image)
