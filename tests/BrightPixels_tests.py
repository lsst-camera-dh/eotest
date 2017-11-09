"""
@brief Unit tests for BrightPixels module.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import unittest
import numpy as np
import lsst.eotest.image_utils as imutils
from lsst.eotest.sensor import MaskedCCD, BrightPixels
import lsst.eotest.sensor.sim_tools as sim_tools


class BrightPixelsTestCase(unittest.TestCase):
    """Test case for BrightPixels code."""

    def setUp(self):
        self.dark_file = 'bright_pixels_test.fits'
        self.emin = 10
        self.exptime = 10
        self.gain = 5
        self.ccdtemp = -100
        self.bias_level, self.bias_sigma = 1e2, 4
        self.dark_curr = 2e-3
        self.ccd = sim_tools.CCD(exptime=self.exptime, gain=self.gain,
                                 ccdtemp=self.ccdtemp)
        self.ccd.add_bias(self.bias_level, self.bias_sigma)
        self.ccd.add_dark_current(level=self.dark_curr)
        self.nsig = ((self.emin*self.exptime/self.gain)
                     / self.ccd.segments[list(self.ccd.segments.keys())[0]].sigma())
        self.npix = 100
        self.pixels = self.ccd.generate_bright_pix(self.npix)
        self.ccd.add_bright_pix(self.pixels, nsig=self.nsig)
        self.ccd.writeto(self.dark_file)

    def tearDown(self):
        os.remove(self.dark_file)

    def test_find_pixel_defects(self):
        ccd = MaskedCCD(self.dark_file)
        for amp in ccd:
            bp = BrightPixels(ccd, amp, ccd.md.get('EXPTIME'),
                              self.gain, ethresh=self.emin/2.,
                              colthresh=100)
            results = bp.find()
            pixels = np.array(np.where(self.pixels[amp] == 1))
            pixels = pixels.transpose()
            pixels = [(x, y) for y, x in pixels]
            pixels.sort()
            self.assertEqual(pixels, results[0])


class BrightColumnsTestCase(unittest.TestCase):
    """Test case for BrightPixels code."""

    def setUp(self):
        self.dark_file = 'bright_columns_test.fits'
        self.emin = 10
        self.exptime = 10
        self.gain = 5
        self.ccdtemp = -100
        self.bias_level, self.bias_sigma = 1e2, 4
        self.dark_curr = 2e-3
        self.ccd = sim_tools.CCD(exptime=self.exptime, gain=self.gain,
                                 ccdtemp=self.ccdtemp)
        self.ccd.add_bias(self.bias_level, self.bias_sigma)
        self.ccd.add_dark_current(level=self.dark_curr)
        self.nsig = ((self.emin*self.exptime/self.gain)
                     / self.ccd.segments[list(self.ccd.segments.keys())[0]].sigma())
        self.ncols = 10
        self.columns = self.ccd.generate_bright_cols(self.ncols)
        self.ccd.add_bright_cols(self.columns, nsig=self.nsig)
        self.ccd.writeto(self.dark_file)

    def tearDown(self):
        os.remove(self.dark_file)

    def test_find_column_defects(self):
        ccd = MaskedCCD(self.dark_file)
        for amp in ccd:
            bp = BrightPixels(ccd, amp, ccd.md.get('EXPTIME'),
                              self.gain, ethresh=self.emin/2.,
                              colthresh=100)
            results = bp.find()
            columns = sorted(self.columns[amp])
            self.assertEqual(columns, results[1])


if __name__ == '__main__':
    unittest.main()
