"""
@brief Find bright pixels and bright columns above a threshold specified
in units of e- per second per pixel.
"""
import os
import numpy as np
import pyfits

import lsst.afw.detection as afwDetect
import lsst.afw.image as afwImage
import lsst.daf.base as dafBase

import lsst.eotest.image_utils as imutils
from MaskedCCD import MaskedCCD
from AmplifierGeometry import makeAmplifierGeometry
from fits_headers import fits_headers

class BrightPixels(object):
    """
    Find bright pixels and bright columns based on a threshold of
    ethresh e- per second per pixel.  A bright column has at least
    colthresh bright pixels.  The bright pixels that are returned are
    exclusive of the bright columns.  The mask that is generated will
    be identified as mask_plane.
    """
    def __init__(self, ccd, amp, exptime, gain,
                 ethresh=5, colthresh=20, mask_plane='BAD'):
        self.ccd = ccd
        self.amp = amp
        self.raw_image = ccd[amp]
        self.exptime = exptime
        self.gain = gain
        self.ethresh = ethresh
        self.colthresh = colthresh
        self.mask_plane = mask_plane
        self.bright_pixels = None
        self.bright_columns = None
    def generate_mask(self, outfile):
        """
        Find bright pixels and columns, and write the mask for this
        amplifier to a FITS file using the afwImage.MaskU.writeFits
        method.
        """
        if self.bright_pixels is None or self.bright_columns is None:
            self.find()
        self.mask = afwImage.MaskU(self.raw_image.getDimensions())
        self.fp_set.setMask(self.mask, self.mask_plane)
        self._write_fits(outfile)
    def _write_fits(self, outfile):
        amp_geom = makeAmplifierGeometry(self.ccd.imfile)
        hdrs = fits_headers()
        if not os.path.isfile(outfile):
            output = pyfits.HDUList()
            output.append(pyfits.PrimaryHDU())
            output[0].header['MASKTYPE'] = 'BRIGHT_PIXELS'
            output[0].header['ETHRESH'] = self.ethresh
            output[0].header['CTHRESH'] = self.colthresh
            output.writeto(outfile, clobber=True)
        imaging = self.ccd.amp_geom.imaging
        ihdr = hdrs[hdrs.keys()[self.amp]]
        md = dafBase.PropertySet()
        md.set('EXTNAME', 'SEGMENT%s' % imutils.channelIds[self.amp])
        md.set('DETSIZE', amp_geom[self.amp]['DETSIZE'])
        md.set('DETSEC', amp_geom[self.amp]['DETSEC'])
        md.set('DATASEC', amp_geom[self.amp]['DATASEC'])
        md.set('NBRTPIX', len(self.bright_pixels))
        md.set('NBRTCOL', len(self.bright_columns))
        self.mask.writeFits(outfile, md, 'a')
    def find(self):
        """
        Find and return the bright pixels and bright columns.
        """
        image = self.ccd.unbiased_and_trimmed_image(self.amp)
        #
        # Multiply e- threshold rate by exptime and convert to DN;
        # create Threshold object.
        #
        threshold = afwDetect.Threshold(self.ethresh*self.exptime/self.gain)
        #
        # Apply footprint detection code.
        #
        self.fp_set = afwDetect.FootprintSet(image, threshold)
        #
        # Organize bright pixels by column.
        #
        columns = dict([(x, []) for x in range(0, self.raw_image.getWidth())])
        for footprint in self.fp_set.getFootprints():
            for span in footprint.getSpans():
                y = span.getY()
                for x in range(span.getX0(), span.getX1()+1):
                    columns[x].append(y)
        #
        # Divide into bright columns (with # bright pixels > self.colthresh)
        # and remaining bright pixels.
        #
        bright_pixs = []
        bright_cols = []
        x0 = image.getX0()
        y0 = image.getY0()
        for x in columns:
            if len(columns[x]) > self.colthresh:
                bright_cols.append(x - x0)
            else:
                bright_pixs.extend([(x - x0, y - y0) for y in columns[x]])
        #
        # Sort the output.
        #
        bright_cols.sort()
        bright_pixs = sorted(bright_pixs)
        self.bright_pixels = bright_pixs
        self.bright_columns = bright_cols
        return bright_pixs, bright_cols
if __name__ == '__main__':
    import sim_tools
    def write_test_image(outfile, emin=10, dark_curr=2e-3, exptime=10,
                         gain=5, ccdtemp=-100, bias_level=1e2,
                         bias_sigma=4, ncols=2, npix=100):
        ccd = sim_tools.CCD(exptime=exptime, gain=gain, ccdtemp=ccdtemp)
        ccd.add_bias(bias_level, bias_sigma)
        ccd.add_dark_current(level=dark_curr)
        #
        # Simulation code sets bright pixel values by nsig*seg.sigma(),
        # but we want to set using e- per sec, so need to compute
        # equivalent nsig.
        #
        nsig = (emin*exptime/gain)/ccd.segments[imutils.allAmps[0]].sigma()
        #
        columns = ccd.generate_bright_cols(ncols)
        ccd.add_bright_cols(columns, nsig=nsig)
        pixels = ccd.generate_bright_pix(npix)
        ccd.add_bright_pix(pixels, nsig=nsig)
        ccd.writeto(outfile)

    def remove_file(filename):
        try:
            os.remove(filename)
        except OSError:
            pass

    dark_file = 'bright_pix_test.fits'
    mask_file = 'bright_pix_mask.fits'
    remove_file(mask_file)

    gain = 5
    write_test_image(dark_file, emin=10, gain=5, npix=1000)

    ccd = MaskedCCD(dark_file)
    for amp in imutils.allAmps:
        bright_pixels = BrightPixels(ccd, amp, ccd.md.get('EXPTIME'), gain)
        pixels, columns = bright_pixels.find()
        bright_pixels.generate_mask(mask_file)
        print imutils.channelIds[amp], len(pixels), columns
