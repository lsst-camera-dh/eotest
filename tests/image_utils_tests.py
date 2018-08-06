"""
@brief Unit tests for image_utils.py module.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import unittest
import numpy as np
import lsst.eotest.image_utils as imutils
from lsst.eotest.sensor import MaskedCCD, AmplifierGeometry
from lsst.eotest.sensor.AmplifierGeometry import makeAmplifierGeometry
import lsst.eotest.sensor.sim_tools as sim_tools
import lsst.afw.image as afwImage

class BiasFunc(object):
    def __init__(self, slope, intercept):
        self.pars = slope, intercept
    def __call__(self, x):
        return x*self.pars[0] + self.pars[1]

class BiasHandlingTestCase(unittest.TestCase):
    bias_slope = 1e-3
    bias_intercept = 0.5
    exptime = 1
    gain = 1
    kwargs = {'fit_order' : 1, 'statistic' : np.mean, 'imaging' : AmplifierGeometry().imaging}
    image_file = 'test_image.fits'
    mean_image_file = 'test_mean_image.fits'
    @classmethod
    def setUpClass(cls):
        cls.amp_geom = AmplifierGeometry()
        cls.bias_image = afwImage.ImageF(cls.amp_geom.full_segment)
        imarr = cls.bias_image.getArray()
        ny, nx = imarr.shape
        yvals = np.arange(0, ny, dtype=np.float)
        bias_func = BiasFunc(cls.bias_slope, cls.bias_intercept)
        for x in range(nx):
            imarr[:,x] += bias_func(yvals)
        ccd = sim_tools.CCD(exptime=cls.exptime, gain=cls.gain,
                            geometry=cls.amp_geom)
        for amp in ccd.segments:
            ccd.segments[amp].image += cls.bias_image
        ccd.writeto(cls.image_file)

        cls.mean_bias_image = afwImage.ImageF(cls.amp_geom.full_segment)
        imarr = cls.mean_bias_image.getArray()
        ny, nx = imarr.shape
        mu, sig = 27316.92, 6.53
        for x in range(nx):
            imarr[:,x] += np.random.normal(mu, sig, ny)
        ccd = sim_tools.CCD(exptime=cls.exptime, gain=cls.gain,
                    geometry=cls.amp_geom)
        for amp in ccd.segments:
            ccd.segments[amp].image += cls.mean_bias_image
        ccd.writeto(cls.mean_image_file)

    @classmethod
    def tearDownClass(cls):
        os.remove(cls.image_file)
    def test_bias_func(self):
        bias_func = BiasFunc(self.bias_slope, self.bias_intercept)
        ccd = MaskedCCD(self.image_file)
        for amp in ccd:
            bf_i = imutils.bias_func(ccd[amp].getImage(),
                                     self.amp_geom.serial_overscan, **self.kwargs)
            bf_m = imutils.bias_func(ccd[amp], self.amp_geom.serial_overscan, **self.kwargs)
            for y in range(2022):
                self.assertEqual(bf_i(y), bf_m(y))
                self.assertAlmostEqual(bias_func(y), bf_m(y))
    def test_bias_image(self):
        ccd = MaskedCCD(self.image_file)
        overscan = makeAmplifierGeometry(self.image_file)
        ccd_mean = MaskedCCD(self.mean_image_file)
        overscan_mean = makeAmplifierGeometry(self.mean_image_file)
        for amp in ccd:
            for method in ['mean', 'row', 'func']:
                if method == 'mean':
		    my_bias_image = imutils.bias_image(ccd_mean[amp],
                                                    overscan_mean.serial_overscan,
                                                    method, **self.kwargs)
                    fracdiff = ((self.mean_bias_image.getArray() - my_bias_image.getArray())
                        /self.mean_bias_image.getArray())
                    self.assertTrue(max(np.abs(fracdiff.flat)) < 1.5e-3)
                else:
		    my_bias_image = imutils.bias_image(ccd[amp],
                                                    overscan.serial_overscan,
                                                    method, **self.kwargs)
                    fracdiff = ((self.bias_image.getArray() - my_bias_image.getArray())
                            /self.bias_image.getArray())
                    self.assertTrue(max(np.abs(fracdiff.flat)) < 1e-6)
    def test_unbias_and_trim(self):
        ccd = MaskedCCD(self.image_file)
        overscan = makeAmplifierGeometry(self.image_file)
        for amp in ccd:
            for method in ['mean', 'row', 'func']:
                image = imutils.unbias_and_trim(ccd[amp], 
                                                 overscan.serial_overscan, method, **self.kwargs)
                imarr = image.getImage().getArray()
                if method == 'mean':
		    self.assertTrue(max(np.abs(imarr.flat)) < 2)
		else:
		    self.assertTrue(max(np.abs(imarr.flat)) < 1e-6)
	    	    #
            	    # Test of corresponding MaskedCCD method.
            	    #
            	    image = ccd.unbiased_and_trimmed_image(amp, overscan.serial_overscan, fit_order=1)
            	    imarr = image.getImage().getArray()
	    	    self.assertTrue(max(np.abs(imarr.flat)) < 1e-6)
    def test_bias_row(self):
        ccd = MaskedCCD(self.image_file)
        overscan = makeAmplifierGeometry(self.image_file)
        for amp in ccd:
            # Unmasked image
            br_i = imutils.bias_row(ccd[amp].getImage(), overscan.serial_overscan)
            # Masked image
            br_m = imutils.bias_row(ccd[amp], overscan.serial_overscan)
	    ny = len(br_i)
	    for ii in range(ny):
	        self.assertEqual(br_i[ii], br_m[ii])
    def test_stack(self):
        ccd = MaskedCCD(self.image_file)
        overscan = makeAmplifierGeometry(self.image_file)
        stats = []
        for method in ['mean', 'row', 'func']:
            corrected = []
            for image in ccd.values():
                corrected.append(imutils.unbias_and_trim(image, overscan.serial_overscan, method, **self.kwargs).getImage())
            stacked = imutils.stack(corrected)
	    imarr = stacked.getArray()
	    if method == 'mean':
		self.assertTrue(max(np.abs(imarr.flat)) < 2)
	    else:
	    	self.assertTrue(max(np.abs(imarr.flat)) < 1e-6)



class FitsMedianTestCase(unittest.TestCase):
    def setUp(self):
        self.values = (0, 1, 2, 3, 4)
        self.files = []
        for i in self.values:
            self.files.append('test_fits_median_image_%02i.fits' % i)
            ccd = sim_tools.CCD(exptime=1)
            for amp in ccd.segments:
                ccd.segments[amp].image += i
            ccd.writeto(self.files[-1])
    def tearDown(self):
        for item in self.files:
            os.remove(item)
    def test_fits_median(self):
        median_image = imutils.fits_median(self.files, hdu=2, fix=True)
        imarr = median_image.getArray()
        for x in imarr.flat:
            self.assertEqual(self.values[len(self.values)/2], x)

if __name__ == '__main__':
    unittest.main()
