"""
@brief Unit tests for MaskedCCD.py module.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import unittest
import numpy as np
import pyfits
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.eotest.image_utils as imutils
import lsst.eotest.sensor as sensorTest
import lsst.eotest.sensor.sim_tools as sim_tools

class BiasFunc(object):
    def __init__(self, slope, intercept):
        self.pars = slope, intercept
    def __call__(self, x):
        return x*self.pars[0] + self.pars[1]

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
        ccd = sim_tools.CCD(exptime=cls.exptime, gain=cls.gain)
        for amp in imutils.allAmps:
            imarr = ccd.segments[amp].image.getArray()
            imarr[cls.ymin:cls.ymax, cls.xmin:cls.xmax] += cls.signal
        ccd.writeto(cls.mask_image)
        cls.mask_files = []
        for mask_plane, bit in cls.mpd.items():
            mask_file = 'mask_file_%s.fits' % mask_plane
            cls.mask_files.append(mask_file)
            masked_ccd = sensorTest.MaskedCCD(cls.mask_image)
            for amp in imutils.allAmps:
                bp = sensorTest.BrightPixels(masked_ccd, amp, cls.exptime,
                                             cls.gain, mask_plane=mask_plane,
                                             ethresh=cls.signal/2.)
                bp.generate_mask(mask_file)
        cls.summed_mask_file = 'summed_mask_file.fits'
        sensorTest.add_mask_files(cls.mask_files, cls.summed_mask_file)
    @classmethod
    def tearDownClass(cls):
        os.remove(cls.mask_image)
        for mask_file in cls.mask_files:
            os.remove(mask_file)
        os.remove(cls.summed_mask_file)
#    @unittest.skip('skip test_add_masks')
    def test_add_masks(self):
        ccd = sensorTest.MaskedCCD(self.mask_image)
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
        ccd = sensorTest.MaskedCCD(self.mask_image)
        for mp, bit in self.mpd.items():
            sctrl = ccd.setMask(mask_name=mp, clear=True)
            self.assertEqual(sctrl.getAndMask(), 2**bit)

class MaskedCCD_biasHandlingTestCase(unittest.TestCase):
    bias_slope = 1e-3
    bias_intercept = 0.5
    exptime = 1
    gain = 1
    image_file = 'test_image.fits'
    @classmethod
    def setUpClass(cls):
        cls.bias_image = afwImage.ImageF(imutils.full_segment)
        imarr = cls.bias_image.getArray()
        ny, nx = imarr.shape
        yvals = np.arange(0, ny, dtype=np.float)
        bias_func = BiasFunc(cls.bias_slope, cls.bias_intercept)
        for x in range(nx):
            imarr[:,x] += bias_func(yvals)
        ccd = sim_tools.CCD(exptime=cls.exptime, gain=cls.gain)
        for amp in imutils.allAmps:
            ccd.segments[amp].image += cls.bias_image
        ccd.writeto(cls.image_file)
    @classmethod
    def tearDownClass(cls):
        os.remove(cls.image_file)
    def test_bias_image(self):
        ccd = sensorTest.MaskedCCD(self.image_file)
        for amp in imutils.allAmps:
            my_bias_image = ccd.bias_image(amp)
            fracdiff = ( (self.bias_image.getArray()-my_bias_image.getArray())
                         /self.bias_image.getArray() )
            self.assertTrue(max(np.abs(fracdiff.flat)) < 1e-6)
    def test_unbias_and_trim(self):
        ccd = sensorTest.MaskedCCD(self.image_file)
        for amp in imutils.allAmps:
            #
            # Test of corresponding MaskedCCD method.
            #
            image = ccd.unbiased_and_trimmed_image(amp)
            imarr = image.getImage().getArray()
            self.assertTrue(max(np.abs(imarr.flat)) < 1e-6)

if __name__ == '__main__':
    unittest.main()
