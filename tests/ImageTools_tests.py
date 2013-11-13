import unittest
import fractions
import numpy as np
import lsst.eotest.utilLib as testUtils
import lsst.afw.image as afwImage

def lcm(*numbers):
    return reduce(lambda x, y : (x*y)/fractions.gcd(x, y), numbers, 1)

class ImageToolsTestCase(unittest.TestCase):
    """Test case for ImageTools static functions"""
    def setUp(self):
        self.binsizes = range(1, 10)
        imsize = lcm(*self.binsizes)
        self.input_image = afwImage.ImageF(imsize, imsize)
        self.input_image += 1
    def tearDow(self):
        pass
    def test_imageTools_rebin(self):
        for binsize in self.binsizes:
            rebinned = testUtils.ImageTools_rebin(self.input_image, binsize)
            rebinned_array = rebinned.getArray()
            self.assertEqual(min(rebinned_array.flat),
                             max(rebinned_array.flat))
            self.assertEqual(min(rebinned_array.flat), binsize**2)

if __name__ == '__main__':
    unittest.main()
