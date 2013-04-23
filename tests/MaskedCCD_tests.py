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
    def setUp(self):
        self.gain = 1
        self.exptime = 1
        self.signal = 1
        self.xmin, self.xmax = 200, 250
        self.ymin, self.ymax = 1000, 1050
        self.mask_image = 'mask_image.fits'
        self.ccd = CCD(self.exptime, gain=self.gain)
        for amp in imutils.allAmps:
            imarr = self.ccd.segments[amp].image.getArray()
            imarr[self.ymin:self.ymax, self.xmin:self.xmax] += self.signal
        self.ccd.writeto(self.mask_image)
        self.mpd = dict(afwImage.MaskU().getMaskPlaneDict().items())
        self.mask_files = []
        for mask_plane, bit in self.mpd.items():
            mask_file = 'mask_file_%s.fits' % mask_plane
            self.mask_files.append(mask_file)
            bp = BrightPixels(self.mask_image, mask_plane=mask_plane,
                              ethresh=self.signal/2.)
            for amp in imutils.allAmps:
                bp.generate_mask(amp, self.gain, mask_file)
        self.summed_mask_file = 'summed_mask_file.fits'
        add_mask_files(self.mask_files, self.summed_mask_file)
    def tearDown(self):
        os.remove(self.mask_image)
        for mask_file in self.mask_files:
            os.remove(mask_file)
        os.remove(self.summed_mask_file)
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
