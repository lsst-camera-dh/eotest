"""
@brief Class to handle masks for a full sensor.  Added mask files can be
or'd to existing masks, and mask bits can be set for use with an
afwMath.makeStatistics object.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import pyfits
import lsst.daf.base as dafBase
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath

import image_utils as imutils

class MaskedCCD(dict):
    def __init__(self, imfile, mask_files=()):
        dict.__init__(self)
        for amp in imutils.allAmps:
            image = afwImage.ImageF(imfile, imutils.dm_hdu(amp))
            mask = afwImage.MaskU(image.getDimensions())
            self[amp] = afwImage.MaskedImageF(image, mask)
        for mask_file in mask_files:
            self.add_masks(mask_file)
        self.stat_ctrl = afwMath.StatisticsControl()
        if mask_files:
            self.setAllMasks()
    def mask_plane_dict(self):
        amp = self.keys()[0]
        return dict(self[amp].getMask().getMaskPlaneDict().items())
    def add_masks(self, mask_file):
        """
        Add a masks from a mask file by or-ing with existing masks.
        """
        for amp in imutils.allAmps:
            curr_mask = self[amp].getMask()
            curr_mask |= afwImage.MaskU(mask_file, imutils.dm_hdu(amp))
    def setMask(self, mask_name=None, clear=False):
        """
        Enable a mask and return the afwMath.StatisticsControl object
        for use by afwMath.makeStatistics. If clear is False, then the
        new mask is or'd with the existing mask.  If clear is False
        and mask_name is None, then all mask bits are cleared.
        """
        if clear:                         # Unset all masks.
            self.stat_ctrl.setAndMask(0)
        if mask_name is not None:         # Set the desired mask.
            new_mask = (self.stat_ctrl.getAndMask()
                        | afwImage.MaskU.getPlaneBitMask(mask_name))
            self.stat_ctrl.setAndMask(new_mask)
        return self.stat_ctrl
    def setAllMasks(self):
        "Enable all masks."
        mpd = self.mask_plane_dict()
        mask_bits = 2**len(mpd) - 1
        self.stat_ctrl.setAndMask(mask_bits)
        return self.stat_ctrl

def add_mask_files(mask_files, outfile, clobber=True):
    masks = dict([(amp, afwImage.MaskU(mask_files[0], imutils.dm_hdu(amp)))
                  for amp in imutils.allAmps])
    for mask_file in mask_files[1:]:
        for amp in imutils.allAmps:
            masks[amp] |= afwImage.MaskU(mask_file, imutils.dm_hdu(amp))
    output = pyfits.HDUList()
    output.append(pyfits.PrimaryHDU())
    output.writeto(outfile, clobber=clobber)
    for amp in imutils.allAmps:
        md = dafBase.PropertySet()
        md.set('EXTNAME', 'AMP%s' % imutils.channelIds[amp])
        md.set('DETSIZE', imutils.detsize)
        md.set('DETSEC', imutils.detsec(amp))
        masks[amp].writeFits(outfile, md, 'a')
    return masks

def compute_stats(image, sctrl, weights=None):
    flags = afwMath.MEAN | afwMath.STDEV
    if weights is None:
        stats = afwMath.makeStatistics(image, flags, sctrl)
    else:
        stats = afwMath.makeStatistics(image, weights, flags, sctrl)
    return stats.getValue(afwMath.MEAN), stats.getValue(afwMath.STDEV)

if __name__ == '__main__':
    image_file = 'bright_pix_test_10.fits'
    mask_files = ('bright_pix_mask.fits', 'CCD250_DEFECTS_mask.fits')

    ccd = MaskedCCD(image_file)
    for mask_file in mask_files:
        print "adding masks from", mask_file
        ccd.add_masks(mask_file)
        print "mask plane dict:", ccd.mask_plane_dict()
        print

    amp = imutils.allAmps[0]
    
    sctrl = ccd.stat_ctrl
    print sctrl.getAndMask(), compute_stats(ccd[amp], sctrl)

    sctrl = ccd.setMask('BAD')
    print sctrl.getAndMask(), compute_stats(ccd[amp], sctrl)

    sctrl = ccd.setMask('CCD250_DEFECTS')
    print sctrl.getAndMask(), compute_stats(ccd[amp], sctrl)

    sctrl = ccd.setMask(clear=True)
    print sctrl.getAndMask(), compute_stats(ccd[amp], sctrl)

    sctrl = ccd.setAllMasks()
    print sctrl.getAndMask(), compute_stats(ccd[amp], sctrl)
