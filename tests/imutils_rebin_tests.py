import unittest
import fractions
import numpy as np
import lsst.eotest.image_utils as imutils
import lsst.afw.image as afwImage


def lcm(*numbers):
    return reduce(lambda x, y: (x*y)/fractions.gcd(x, y), numbers, 1)


class RebinTestCase(unittest.TestCase):
    """Test case for image_utils.rebin function."""

    def setUp(self):
        self.binsizes = range(1, 10)
        imsize = lcm(*self.binsizes)
        self.input_image = afwImage.ImageF(imsize, imsize)
        self.input_image += 1

    def tearDown(self):
        pass

    def test_image_utils_rebin(self):
        for binsize in self.binsizes:
            rebinned = imutils.rebin(self.input_image, binsize)
            rebinned_array = rebinned.getArray()
            self.assertEqual(min(rebinned_array.flat),
                             max(rebinned_array.flat))
            self.assertEqual(min(rebinned_array.flat), binsize**2)


if __name__ == '__main__':
    unittest.main()
