import numpy as np
import pyfits

import lsst.afw.display.ds9 as ds9
import lsst.afw.detection as afwDetect
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath

import image_utils as imutils

class MaskedCCD(dict):
    def __init__(self, imfile):
        dict.__init__(self)
        for amp in imutils.allAmps:
            image = afwImage.ImageF(imfile, imutils.dm_hdu(amp))
            mask = afwImage.MaskU(image.getDimensions())
            self[amp] = afwImage.MaskedImageF(image, mask)
    def add_mask(self, mask_file, mask_name='BAD', threshold=0.5):
        mpd = self._update_mask_planes(mask_name,
                                       mask=self[imutils.allAmps[0]].getMask())
        threshold = afwDetect.Threshold(threshold)
        for amp in imutils.allAmps:
            my_mask = afwImage.ImageU(mask_file, imutils.dm_hdu(amp))
            #
            # MaskU.setMaskPlaneValues seems to want Spans from
            # afwDetect.Footprints for setting the mask plane values.
            #
            fp_set = afwDetect.FootprintSet(my_mask, threshold)
            for fp in fp_set.getFootprints():
                for span in fp.getSpans():
                    self[amp].getMask().setMaskPlaneValues(mpd[mask_name],
                                                           span.getX0(),
                                                           span.getX1(),
                                                           span.getY())
    def _update_mask_planes(self, mask_name, mask=None):
        if mask is None:
            mask = afwImage.MaskU(1, 1)
        mpd = dict(mask.getMaskPlaneDict().items())
        if not mpd.has_key(mask_name):
            mask.addMaskPlane(mask_name)
        return mpd
    def setAndMask(self, mask_name=None, clear=False):
        stat_ctrl = afwMath.StatisticsControl()
        if clear:     # Unset all masks.
            stat_ctrl.setAndMask(0)
        if mask_name is not None: # Set the desired mask.
            stat_ctrl.setAndMask(afwImage.MaskU.getPlaneBitMask(mask_name))
        return stat_ctrl

if __name__ == '__main__':           
    image_file = '001-00/dark/data/001_00_dark_100_000.fits'
    mask_file = '001-00/bright_pixels/data/001-00_bright_pixel_map.fits'
    mask_name = 'BAD'

    flags = afwMath.MEAN | afwMath.STDEV

    ccd = MaskedCCD(image_file)
    ccd.add_mask(mask_file, mask_name)

    #
    # Compute statistics with bad pixels masked out.
    #
    stat_ctrl = ccd.setAndMask(mask_name)
    for amp in imutils.allAmps:
        stats = afwMath.makeStatistics(ccd[amp], flags, stat_ctrl)
        print amp, stats.getValue(afwMath.MEAN), stats.getValue(afwMath.STDEV)
    print
    
    #
    # Unset bad-pixel mask and recompute statistics
    #
    stat_ctrl = ccd.setAndMask(clear=True)
    for amp in imutils.allAmps:
        stats = afwMath.makeStatistics(ccd[amp], flags, stat_ctrl)
        print amp, stats.getValue(afwMath.MEAN), stats.getValue(afwMath.STDEV)
