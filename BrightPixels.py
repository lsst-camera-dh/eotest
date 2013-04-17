"""
@brief Find bright pixels and bright columns above a threshold specified
in units of e- per second per pixel.
"""
import os
import pyfits
import lsst.afw.detection as afwDetect
import lsst.afw.image as afwImage
import lsst.daf.base as dafBase
import image_utils as imutils
from simulation.sim_tools import CCD

class BrightPixels(object):
    """
    Find bright pixels and bright columns based on a threshold of 5e-
    per second per pixel.  A bright column has at least 20 bright
    pixels.  The bright pixels that are returned are exclusive of the
    bright columns.
    """
    def __init__(self, dark_file, amp):
        self.dark_file = dark_file
        self.amp = amp
        self.md = afwImage.readMetadata(dark_file, 1)
        self.raw_image = afwImage.ImageF(dark_file, imutils.dm_hdu(amp))
        self.mask = afwImage.MaskU(self.raw_image.getDimensions())
    def write_mask(self, outfile, mask_plane='BAD'):
        """
        Write the mask for this amplifier to a FITS file using the
        afwImage.MaskU.writeFits method.  This method writes the mask
        plane bit mapping as header keywords like MP_BAD = 0, MP_SAT =
        1, etc.. Also add keywords for iraf mosaicking in ds9.
        """
        self.fp_set.setMask(self.mask, mask_plane)
        if not os.path.isfile(outfile):
            foo = pyfits.HDUList()
            foo.append(pyfits.PrimaryHDU())
            foo.writeto(outfile, clobber=True)
        md = dafBase.PropertySet()
        md.set('EXTNAME', 'AMP%s' % imutils.channelIds[self.amp])
        md.set('DETSIZE', imutils.detsize)
        md.set('DETSEC', imutils.detsec(self.amp))
        self.mask.writeFits(outfile, md, 'a')
    def find(self, gain, ethresh=5, colthresh=20):
        """
        Find and return the bright pixels and bright columns.
        """
        exptime = self.md.get('EXPTIME')
        #
        image = imutils.unbias_and_trim(self.raw_image)
        #
        # Multiply e- threshold rate by exptime and convert to DN;
        # create Threshold object.
        #
        threshold = afwDetect.Threshold(ethresh*exptime/gain)
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
        # Divide into bright columns (with # bright pixels > colthresh) and
        # remaining bright pixels.
        #
        bright_pixs = []
        bright_cols = []
        for y in columns:
            if len(columns[y]) > colthresh:
                bright_cols.append(y)
            else:
                bright_pixs.extend([(x, y) for x in columns[y]])
        #
        # Sort the output.
        #
        bright_cols.sort()
        bright_pixs = sorted(bright_pixs)
        return bright_pixs, bright_cols

def write_test_image(outfile, emin=10, dark_curr=2e-3, exptime=10,
                     gain=5, ccdtemp=-100, bias_level=1e2,
                     bias_sigma=4, ncols=2, npix=100):
    ccd = CCD(exptime=exptime, gain=gain, ccdtemp=ccdtemp)
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

if __name__ == '__main__':
    dark_file = 'bright_pix_test_10.fits'
    mask_file = 'bright_pix_mask.fits'
    remove_file(mask_file)

    gain = 5
    #write_test_image(dark_file, emin=10, gain=5, npix=1000)
    
    for amp in imutils.allAmps:
        bright_pixels = BrightPixels(dark_file, amp)
        pixels, columns = bright_pixels.find(gain, ethresh=5)
        print imutils.channelIds[amp], len(pixels), columns
        bright_pixels.write_mask(mask_file)
