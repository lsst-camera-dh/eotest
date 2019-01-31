"""
@brief Unit tests for BFTask module.
"""
import os
import unittest
import numpy as np
import lsst.eotest.image_utils as imutils
from lsst.eotest.sensor import MaskedCCD
import lsst.eotest.sensor.sim_tools as sim_tools
from lsst.eotest.sensor.BFTask import crossCorrelate


class BFTestCase(unittest.TestCase):
    """Test case for BrightFatter code.
        Modeled after BrightPixelsTestCase"""

    def setUp(self):
        self.flat_file1 = 'brighter_fatter_test1.fits'
        self.flat_file2 = 'brighter_fatter_test2.fits'
        self.emin = 10
        self.exptime = 10
        self.gain = 5
        self.ccdtemp = -100
        self.maxLag = 1
        self.nPixBorder = 10
        self.nSigmaClip = 3
        self.backgroundBinSize = 128
        self.bias_level, self.bias_sigma = 1e2, 4
        self.ccd1 = sim_tools.CCD(exptime=self.exptime, gain=self.gain,
                                 ccdtemp=self.ccdtemp)
        self.ccd1.add_bias(self.bias_level, self.bias_sigma)

        self.ccd1.writeto(self.flat_file1)
        self.ccd2 = sim_tools.CCD(exptime=self.exptime, gain=self.gain,
                                 ccdtemp=self.ccdtemp)
        self.ccd2.add_bias(self.bias_level, self.bias_sigma)

        self.ccd2.writeto(self.flat_file2)

    def tearDown(self):
        os.remove(self.flat_file1)
        os.remove(self.flat_file2)

    def test_find_pixel_defects(self):
        ccd1 = MaskedCCD(self.flat_file1)
        ccd2 = MaskedCCD(self.flat_file2)
        for amp in ccd1:
            xcorr, xcorr_err = crossCorrelate(ccd1.unbiased_and_trimmed_image(amp), ccd2.unbiased_and_trimmed_image(amp), self.maxLag, self.nSigmaClip, self.backgroundBinSize )
            self.assertTrue(1>xcorr[0][1] and 1>xcorr[1][0])


if __name__ == '__main__':
    unittest.main()
