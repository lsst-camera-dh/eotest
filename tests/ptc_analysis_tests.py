"""
Test code for PTC analysis.
"""
import os
import unittest
import numpy as np
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
        """
        Run the fit using a canned file representing mean-var for a
        single amplifier.
        """
        infile = os.path.join(os.environ['EOTEST_DIR'], 'tests',
                              'PTC_mean_var_values.txt')
        data = np.loadtxt(infile)
        mean = data[:, 0]
        var = data[:, 1]

        # Fit mean-variance relation, compare to expected results
        # for gain, noise, k, rolloff.
        (ptc_gain, ptc_gain_error, ptc_a00, ptc_a00_error, ptc_noise,
         ptc_noise_error, ptc_turnoff) \
         = sensorTest.PtcTask.fit_ptc_curve(mean, var)

        self.assertAlmostEqual(ptc_gain, 0.77132744, places=5)
        self.assertAlmostEqual(ptc_gain_error, 0.001864, places=5)
        self.assertAlmostEqual(ptc_a00, 2.67726e-6, places=5)
        self.assertAlmostEqual(ptc_a00_error, 3.94e-8, places=3)
        self.assertAlmostEqual(ptc_noise, 5.2091694, places=4)
        self.assertAlmostEqual(ptc_noise_error, 0.38651699, places=4)
        self.assertAlmostEqual(ptc_turnoff, 148175, places=3)

if __name__ == '__main__':
    unittest.main()
