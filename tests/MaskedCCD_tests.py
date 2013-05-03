"""
@brief Unit tests for MaskedCCD.py module.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import unittest
import pyfits
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import image_utils as imutils
from simulation.sim_tools import CCD
from BrightPixels import BrightPixels
from MaskedCCD import MaskedCCD, add_mask_files

class MaskedCCDTestCase(unittest.TestCase):
    gain = 1
    exptime = 1
    signal = 1
    xmin, xmax = 200, 250
    ymin, ymax = 1000, 1050
    mask_image = 'mask_image.fits'
    mpd = dict(afwImage.MaskU().getMaskPlaneDict().items())
    @classmethod
    def setUpClass(cls):
        ccd = CCD(exptime=cls.exptime, gain=cls.gain)
        for amp in imutils.allAmps:
            imarr = ccd.segments[amp].image.getArray()
            imarr[cls.ymin:cls.ymax, cls.xmin:cls.xmax] += cls.signal
        ccd.writeto(cls.mask_image)
        cls.mask_files = []
        for mask_plane, bit in cls.mpd.items():
            mask_file = 'mask_file_%s.fits' % mask_plane
            cls.mask_files.append(mask_file)
            masked_ccd = MaskedCCD(cls.mask_image)
            for amp in imutils.allAmps:
                bp = BrightPixels(masked_ccd[amp], cls.exptime, cls.gain,
                                  mask_plane=mask_plane, ethresh=cls.signal/2.)
                bp.generate_mask(mask_file, amp)
        cls.summed_mask_file = 'summed_mask_file.fits'
        add_mask_files(cls.mask_files, cls.summed_mask_file)
    @classmethod
    def tearDownClass(cls):
        os.remove(cls.mask_image)
        for mask_file in cls.mask_files:
            os.remove(mask_file)
        os.remove(cls.summed_mask_file)
#    @unittest.skip('skip test_add_masks')
    def test_add_masks(self):
        ccd = MaskedCCD(self.mask_image)
        ccd.add_masks(self.summed_mask_file)
        total_signal = (self.signal*(self.ymax - self.ymin)
                        *(self.xmax - self.xmin))
        ny, nx = ccd[imutils.allAmps[0]].getImage().getArray().shape
        for amp in imutils.allAmps:
            self.assertEqual(sum(ccd[amp].getImage().getArray().flat),
                             total_signal)
        stat_ctrl = ccd.setMask('BAD')
        for amp in imutils.allAmps:
            stats = afwMath.makeStatistics(ccd[amp], afwMath.MEAN, stat_ctrl)
            self.assertEqual(0, stats.getValue(afwMath.MEAN))
        stat_ctrl = ccd.setMask(clear=True)
        for amp in imutils.allAmps:
            stats = afwMath.makeStatistics(ccd[amp], afwMath.MEAN, stat_ctrl)
            self.assertAlmostEqual(float(total_signal)/(nx*ny),
                                   stats.getValue(afwMath.MEAN), places=10)
#    @unittest.skip('skip test_setMask')
    def test_setMask(self):
        ccd = MaskedCCD(self.mask_image)
        for mp, bit in self.mpd.items():
            sctrl = ccd.setMask(mask_name=mp, clear=True)
            self.assertEqual(sctrl.getAndMask(), 2**bit)

if __name__ == '__main__':
    unittest.main()
