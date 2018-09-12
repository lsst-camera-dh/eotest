"""
Test cte_matrix using EPER calculation.
"""
import unittest
import numpy as np
import lsst.eotest.sensor as sensorTest


class Cte_Tools_TestCase(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_cte_matrix(self):
        self.cti_values = np.logspace(-7, -4, 16)
        signal = 1e4
        npix = 2000
        noscan = 20
        for i, cti in enumerate(self.cti_values):
            pixels = np.zeros(npix, dtype=np.float)
            pixels[:-noscan] += signal
            cte = sensorTest.cte_matrix(len(pixels), cti)
            pixels = np.dot(cte, pixels)
            cti_est = \
                sum(pixels[npix-noscan:])/pixels[npix-noscan-1]/(npix-noscan+1)
            self.assertAlmostEqual(1, cti_est/cti)


if __name__ == '__main__':
    unittest.main()
