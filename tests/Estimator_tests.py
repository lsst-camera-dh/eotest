"""
Unit tests for Estimator class.
"""
import unittest
import numpy as np
import numpy.random as ra
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
from lsst.eotest.Estimator import Estimator


class EstimatorTestCase(unittest.TestCase):
    def setUp(self):
        np.seterr('raise')
        ra.seed(930293401)
        nx, ny = 100, 100
        self.mean1, self.sigma1 = 1, 0.1
        self.mean2, self.sigma2 = 2, 0.3
        # The Estimator class should work for both Image and MaskedImage
        self.image1 = afwImage.ImageF(nx, ny)
        self.image1.getArray()[:] = ra.normal(self.mean1, self.sigma1, (nx, ny))
        self.image2 = afwImage.MaskedImageF(nx, ny)
        self.image2.getImage().getArray()[:] = ra.normal(self.mean2,
                                                         self.sigma2, (nx, ny))
        self.stat_ctrl = None
        self.gain = 2
        self.est1 = Estimator(self.image1, self.stat_ctrl, gain=self.gain)
        self.est2 = Estimator(self.image2, self.stat_ctrl, gain=self.gain,
                              statistic=afwMath.MEDIAN)

    def tearDown(self):
        pass

    def test_stats(self):
        arr1 = self.image1.getArray()*self.gain
        value1 = np.mean(arr1)
        error1 = np.sqrt(sum(arr1.flat))/np.prod(arr1.shape)
        self.assertAlmostEqual(value1, self.est1.value, places=6)
        self.assertAlmostEqual(error1, self.est1.error, places=6)

        arr2 = self.image2.getImage().getArray()*self.gain
        value2 = np.median(arr2)
        error2 = np.sqrt(sum(arr2.flat))/np.prod(arr2.shape)
        self.assertAlmostEqual(value2, self.est2.value, places=6)
        self.assertAlmostEqual(error2, self.est2.error, places=6)

    def test_addition(self):
        result = self.est1 + self.est2
        self.assertEqual(result.value, self.est1.value + self.est2.value)
        self.assertEqual(result.error, np.sqrt(self.est1.error**2 +
                                               self.est2.error**2))
        foo = sum([self.est1, self.est2])
        self.assertEqual(foo.value, result.value)
        self.assertEqual(foo.error, result.error)

    def test_subtraction(self):
        result = self.est1 - self.est2
        self.assertEqual(result.value, self.est1.value - self.est2.value)
        self.assertEqual(result.error, np.sqrt(self.est1.error**2 +
                                               self.est2.error**2))

    def test_rsubtraction(self):
        my_term = 10.
        result = my_term - self.est1
        self.assertEqual(result.value, my_term - self.est1.value)
        self.assertEqual(result.error, self.est1.error)

    def test_multiplication(self):
        my_value = self.est1.value*self.est2.value
        my_error = np.abs(my_value) * \
            np.sqrt((self.est1.error/self.est1.value)**2 +
                    (self.est2.error/self.est2.value)**2)
        result = self.est1*self.est2
        self.assertEqual(result.value, my_value)
        self.assertEqual(result.error, my_error)

    def test_division(self):
        my_value = self.est1.value/self.est2.value
        my_error = np.abs(my_value) * \
            np.sqrt((self.est1.error/self.est1.value)**2 +
                    (self.est2.error/self.est2.value)**2)
        result = self.est1/self.est2
        self.assertEqual(result.value, my_value)
        self.assertEqual(result.error, my_error)

    def test_equation(self):
        const = 5
        result = (self.est1 + const*self.est2)/self.est2
        numerator_value = self.est1.value + const*self.est2.value
        numerator_err = np.sqrt(self.est1.error**2 + (const*self.est2.error)**2)
        denominator_value = self.est2.value
        denominator_err = self.est2.error
        my_value = numerator_value/denominator_value
        my_error = np.abs(my_value) * \
            np.sqrt((numerator_err/numerator_value)**2 +
                    (denominator_err/denominator_value)**2)
        self.assertAlmostEqual(result.value, my_value)
        self.assertAlmostEqual(result.error, my_error)


if __name__ == '__main__':
    unittest.main()
