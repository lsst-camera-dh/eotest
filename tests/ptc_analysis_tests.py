"""
Test code for PTC analysis.
"""
import os
import sys
import unittest
import numpy as np
import astropy.io.fits as fits
import lsst.eotest.sensor as sensorTest

class PTCGainFitterTestCase(unittest.TestCase):
    """
    TestCase class for PTC gain fitting, plus k factor, noise, and turn-over
    """
    def setUp(self):
        # Save current floating point error handling.
        self.np_fp_config = np.geterr()
        # Explicitly ignore underflows.
        np.seterr(under='ignore')

    def tearDown(self):
        # Restore floating point error handling.
        np.seterr(**self.np_fp_config)

    def test_ptc_fit(self):
        # Run the fit using a canned file representing mean-var for a
        # single amplifier.
        infile = os.path.join(os.environ['EOTEST_DIR'], 'tests',
                              'PTC_mean_var_values.txt')
        data = np.loadtxt(infile)
        mean = data[:, 0]
        var = data[:, 1]
        ptc_stats = {}
        ptc_stats[0] = [mean, var]

        # Fit mean-variance relation, compare to expected results
        # for gain, noise, k, rolloff
        # Below writes a file TEST_ID_eotest_results.fits
        task = sensorTest.PtcTask()
        task._fit_curves(ptc_stats, 'TEST_ID')
        eotest_file = 'TEST_ID_eotest_results.fits'

        # Read the file and check the results
        with fits.open(eotest_file) as hdulist:
            hdudata = hdulist[1].data
            self.assertAlmostEqual(hdudata['PTC_GAIN'][0], 0.77132744, places=5)
            self.assertAlmostEqual(hdudata['PTC_GAIN_ERROR'][0], 0.001864, places=5)
            self.assertAlmostEqual(hdudata['PTC_A00'][0], 2.67726e-6, places=5)
            self.assertAlmostEqual(hdudata['PTC_A00_ERROR'][0], 3.94e-8, places=3)
            self.assertAlmostEqual(hdudata['PTC_NOISE'][0], 5.2091694, places=4)
            self.assertAlmostEqual(hdudata['PTC_NOISE_ERROR'][0], 0.38651699, places=4)
            self.assertAlmostEqual(hdudata['PTC_TURNOFF'][0], 148175, places=3)

        os.remove(eotest_file)

if __name__ == '__main__':
    unittest.main()
