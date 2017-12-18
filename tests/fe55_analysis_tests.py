"""
Test code for Fe55 analyses.
"""
from __future__ import print_function
import os
import unittest
import numpy as np
import lsst.eotest.sensor as sensorTest

class Fe55GainFitterTestCase(unittest.TestCase):
    def setUp(self):
        # Save current floating point error handling.
        self.np_fp_config = np.geterr()
        # Explicitly ignore underflows.
        np.seterr(under='ignore')

    def tearDown(self):
        # Restore floating point error handling.
        np.seterr(**self.np_fp_config)

    def test_set_hist_range(self):
        infile = os.path.join(os.environ['EOTEST_DIR'], 'tests',
                              'Fe55_DN_values.txt')
        dn0 = np.array(np.recfromtxt(infile))
        mode0 = 2300.
        gain0 = 0.69
        dADU = 50
        bins = 100
        hist_nsig = 10
        for gain in np.logspace(-1, 1, 20):
            mode = gain0/gain*mode0
            dn = gain0/gain*dn0
            fitter = sensorTest.Fe55GainFitter(dn)
            fitter._set_hist_range(dADU, bins, hist_nsig)
            self.assertLess(fitter.xrange[0], mode)
            self.assertLess(mode, fitter.xrange[1])
            fitter.fit()

if __name__ == '__main__':
    unittest.main()
