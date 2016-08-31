"""
@brief Unit tests for DarkPixels module.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import unittest
import numpy as np
import lsst.eotest.image_utils as imutils
from lsst.eotest.sensor import MaskedCCD, DarkPixels
import lsst.eotest.sensor.sim_tools as sim_tools

class DarkPixelsTestCase(unittest.TestCase):
    """Test case for DarkPixels code."""
    def setUp(self):
        self.sflat_file = 'dark_pixels_test.fits'
        self.emin = 10
        self.exptime = 10
        self.gain = 5
        self.ccdtemp = -100
        self.bias_level = 1e2
        self.bias_sigma = 4
        self.dark_curr = 2e-3
        self.frac_level = 0.5
        self.ccd = sim_tools.CCD(exptime=self.exptime, gain=self.gain,
                                 ccdtemp=self.ccdtemp)
        self.ccd.add_bias(self.bias_level, self.bias_sigma)
        self.ccd.add_dark_current(level=self.dark_curr)
        # For an exptime of 10, this should give about 2e5 e- per pixel.
        self.ccd.expose_flat(1e4)
        # Number of dark pixels per segment
        self.npix = 100
        # This just generates a mask of pixel locations for each amp.
        self.pixels = self.ccd.generate_bright_pix(self.npix)
        self.ccd.set_dark_pix(self.pixels, self.frac_level)
        self.ccd.writeto(self.sflat_file)
    def tearDown(self):
        os.remove(self.sflat_file)
    def test_find_pixel_defects(self):
        ccd = MaskedCCD(self.sflat_file)
        for amp in ccd:
            dp = DarkPixels(ccd, amp, frac_thresh=0.8, colthresh=100)
            results = dp.find()
            pixels = np.array(np.where(self.pixels[amp] == 1))
            pixels = pixels.transpose()
            pixels = [(x, y) for y, x in pixels]
            pixels.sort()
            self.assertEqual(pixels, results[0])

class DarkColumnsTestCase(unittest.TestCase):
    """Test case for DarkPixels code."""
    def setUp(self):
        self.sflat_file = 'dark_columns_test.fits'
        self.emin = 10
        self.exptime = 10
        self.gain = 5
        self.ccdtemp = -100
        self.bias_level = 1e2
        self.bias_sigma = 4
        self.dark_curr = 2e-3
        self.frac_level = 0.5
        self.ccd = sim_tools.CCD(exptime=self.exptime, gain=self.gain,
                                 ccdtemp=self.ccdtemp)
        self.ccd.add_bias(self.bias_level, self.bias_sigma)
        self.ccd.add_dark_current(level=self.dark_curr)
        # For an exptime of 10, this should give about 2e5 e- per pixel.
        self.ccd.expose_flat(1e4)
        self.ncols = 10
        self.columns = self.ccd.generate_bright_cols(self.ncols)
        self.ccd.set_dark_cols(self.columns, self.frac_level)
        self.ccd.writeto(self.sflat_file)
    def tearDown(self):
        os.remove(self.sflat_file)
    def test_find_column_defects(self):
        ccd = MaskedCCD(self.sflat_file)
        for amp in ccd:
            dp = DarkPixels(ccd, amp, frac_thresh=0.8, colthresh=100)
            results = dp.find()
            columns = sorted(self.columns[amp])
            self.assertEqual(len(columns), self.ncols)
            self.assertEqual(len(results[1]), self.ncols)
            self.assertEqual(columns, results[1])

if __name__ == '__main__':
    unittest.main()
