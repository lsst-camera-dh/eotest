"""
Test code for PTC analysis.
"""
from __future__ import print_function
import os
import unittest
import itertools
import numpy as np
import lsst.eotest.sensor as sensorTest
import astropy.io.fits as fits
import sys

class PTCGainFitterTestCase(unittest.TestCase):
    def setUp(self):
        # Save current floating point error handling.
        self.np_fp_config = np.geterr()
        # Explicitly ignore underflows.
        np.seterr(under='ignore')

    def tearDown(self):
        # Restore floating point error handling.
        np.seterr(**self.np_fp_config)

    def test_ptc_fit(self):
        infile = os.path.join(os.environ['EOTEST_DIR'], 'tests',
                              'PTC_mean_var_values.txt')
        data = np.loadtxt(infile)
        mean = data[0]
        var = data[1]
        ptc_stats = {}
        ptc_stats['Amp0'] = [mean,var]

        # Fit mean-variance relation, compare to expected results
        # for gain, noise, k, rolloff
        # Below writes a file TEST_ID_eotest_results.fits
        ptc_task._fit_curves(ptc_stats, 'TEST_ID')

        # Read the file and check the results
        try:
            hdu = fits.open('TEST_ID_eotest_results.fits')
        except:
            print('File not found:  TEST_ID_eotest_results.fits')
            sys.exit()

        self.assertAlmostEqual(hdu[1].data['PTC_GAIN'][0], 0.767, places=3)
        self.assertAlmostEqual(hdu[1].data['PTC_GAIN_ERROR'][0], 0.002, places=3)
        self.assertAlmostEqual(hdu[1].data['PTC_A00'][0], 2.68e-6, places=3)
        self.assertAlmostEqual(hdu[1].data['PTC_A00_ERROR'][0], 3.94e-8, places=3)
        self.assertAlmostEqual(hdu[1].data['PTC_NOISE'][0], 5.0, places=1)
        self.assertAlmostEqual(hdu[1].data['PTC_NOISE_ERROR'][0], 0.4, places=1)
        self.assertAlmostEqual(hdu[1].data['PTC_TURNOFF'][0], 150000, places=3)

if __name__ == '__main__':
    unittest.main()
