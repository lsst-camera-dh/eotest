"""
Test code for Fe55 analyses.
"""
from __future__ import print_function
import os
import unittest
import itertools
import numpy as np
import lsst.eotest.sensor as sensorTest

class Fe55GainFitterTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

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
            #fitter.fit()

class Fe55PsfTestCase(unittest.TestCase):
    "Test case class for fe55_psf module"
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_p9_values(self):
        "Test evaluation of the p9 pixel values."
        class Peak(object):
            "Proxy class for footprint peak object."
            def __init__(self, ix, iy):
                "Constructor"
                self.ix = ix
                self.iy = iy
            def getIx(self):
                "Return x-pixel index."
                return self.ix
            def getIy(self):
                "Return y-pixel index."
                return self.iy
        peak = Peak(1, 1)
        imarr = np.arange(9).reshape((3, 3))
        x0, y0 = 1, 1
        sigmax, sigmay = 0.5, 0.5
        DN_tot = 9.
        pos = [(x[1], x[0]) for x in itertools.product((-1, 0, 1), (-1, 0, 1))]
        p9_data, p9_model \
            = sensorTest.fe55_psf.p9_values(peak, imarr, x0, y0,
                                            sigmax, sigmay, DN_tot)
        self.assertEqual(len(p9_data), len(p9_model))
        self.assertItemsEqual(p9_data, range(9))
        for i in (2, 6, 8):
            self.assertEqual(p9_model[0], p9_model[i])
        for i in (3, 5, 7):
            self.assertEqual(p9_model[1], p9_model[i])
        for i in range(4) + range(5, 9):
            self.assertGreater(p9_model[4], p9_model[i])

if __name__ == '__main__':
    unittest.main()
