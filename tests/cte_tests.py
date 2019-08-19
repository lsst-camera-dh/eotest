"""
Test cte_matrix using EPER calculation.
"""
import unittest
import numpy as np
import scipy.special
import lsst.eotest.sensor as sensorTest


class Cte_Tools_TestCase(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_cte_matrix(self):
        """Test that the inferred CTI using EPER agrees with input value."""
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
                sum(pixels[-noscan:])/pixels[npix-noscan-1]/(npix-noscan)
            cti_est /= (1 + cti_est)
            #print(cti, cti_est, cti_est/cti)
            self.assertAlmostEqual(1, cti_est/cti)

    def test_3pixel_matrix(self):
        """Test against a 3-pixel CTE matrix derived by hand."""
        npix = 3
        cti = 1e-5
        cte = sensorTest.cte_matrix(npix, cti)
        # Test individual matrix element values.
        cte_ref = np.array([(         1-cti,                0,          0),
                            (   cti*(1-cti),       (1-cti)**2,          0),
                            (cti**2*(1-cti), 2*cti*(1-cti)**2, (1-cti)**3)])

        for _ in zip(cte.ravel(), cte_ref.ravel()):
            self.assertAlmostEqual(*_)

if __name__ == '__main__':
    unittest.main()
