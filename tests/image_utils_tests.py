import os
import unittest
import numpy as np
import lsst.afw.image as afwImage
from simulation.sim_tools import CCD
from MaskedCCD import MaskedCCD
import image_utils as imutils

class BiasHandlingTestCase(unittest.TestCase):
    def setUp(self):
        self.bias_slope = 1e-3
        self.bias_intercept = 0.5
        self.exptime = 1
        self.gain = 1
        self.image_file = 'test_image.fits'

        self.bias_func = lambda pixel : pixel*self.bias_slope \
                         + self.bias_intercept
        self.bias_image = afwImage.ImageF(imutils.full_segment)
        imarr = self.bias_image.getArray()
        ny, nx = imarr.shape
        yvals = np.arange(0, ny, dtype=np.float)
        for x in range(nx):
            imarr[:,x] += self.bias_func(yvals)

        ccd = CCD(exptime=self.exptime, gain=self.gain)
        for amp in imutils.allAmps:
            ccd.segments[amp].image += self.bias_image
        ccd.writeto(self.image_file)
    def tearDown(self):
        os.remove(self.image_file)
    def test_bias_func(self):
        ccd = MaskedCCD(self.image_file)
        for amp in imutils.allAmps:
            bf_i = imutils.bias_func(ccd[amp].getImage())
            bf_m = imutils.bias_func(ccd[amp])
            for y in range(2022):
                self.assertEqual(bf_i(y), bf_m(y))
                self.assertAlmostEqual(self.bias_func(y), bf_m(y))
    def test_bias_image(self):
        ccd = MaskedCCD(self.image_file)
        for amp in imutils.allAmps:
            my_bias_image = imutils.bias_image(ccd[amp])
            fracdiff = ( (self.bias_image.getArray()-my_bias_image.getArray())
                         /self.bias_image.getArray() )
            self.assertTrue(max(np.abs(fracdiff.flat)) < 1e-6)
    def test_unbias_and_trim(self):
        ccd = MaskedCCD(self.image_file)
        for amp in imutils.allAmps:
            image = imutils.unbias_and_trim(ccd[amp])
            imarr = image.getImage().getArray()
            self.assertTrue(max(np.abs(imarr.flat)) < 1e-6)

if __name__ == '__main__':
    unittest.main()

