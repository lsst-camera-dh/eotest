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
    def __init__(self, exptime=1, gain=5, bbox=full_segment):
        self.exptime = exptime
        self.gain = gain
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
    def add_dark_current(self, level=10):
        dark_arr = self._poisson_imarr(level*self.exptime)
        self.imarr += dark_arr
    def expose_flat(self, level):
        flat_arr = self._poisson_imarr(level*self.exptime)*self.gain
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
    def add_Fe55_hits(self, nxrays=200):
        """One- and two-pixel distributions, based on statistics
        from one of PD's Fe55 exposures."""
        ny, nx = self.imarr.shape
        for i in range(nxrays):
            x0 = random.randint(nx)
            y0 = random.randint(ny)
            signal = random.normal(1620., 13.)/self.gain  # Janesick, p 132.
            if random.uniform() < 0.5:  # Single pixel hit
                self.imarr[y0][x0] += int(signal)
            else:  # Two pixels
                peak_ratio = min(random.normal(0.58, 0.055), 1)
                self.imarr[y0][x0] = int(signal*peak_ratio)
                tail = int(1 - signal*peak_ratio)
                try:
                    if random.uniform() < 0.5:
                        self.imarr[y0][x0 + random.randint(2)*2 - 1] += tail
                    else:
                        self.imarr[y0 + random.randint(2)*2 - 1][x0] += tail
                except IndexError:
                    pass
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
    return outfile
    
def simulateDark(outfile, dark_curr, exptime=1, hdus=16, verbose=True):
    if verbose:
        print "simulating dark:", outfile
    segments = []
    for i in range(hdus):
        if verbose:
            print "HDU", i
        seg = SegmentExposure(exptime=exptime)
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
        seg = SegmentExposure(exptime=exptime, gain=gain)
        seg.add_bias()
        seg.add_dark_current(dark_curr)
        seg.expose_flat(level)
        segments.append(seg)
    writeFits(segments, outfile)

if __name__ == '__main__':
    seg = SegmentExposure()
    
    seg.add_bias(1e4, 10)
    seg.add_dark_current(300)
    seg.expose_flat(200)
    cols = seg.add_bright_cols(ncols=1, nsig=5)
    pix = seg.add_bright_pix(npix=100, nsig=10)

    writeFits((seg,), 'test_image.fits')

#    ds9.mtv(exp.image)
