"""
@brief Unit tests for MaskedCCD.py module.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import unittest
import numpy as np
import lsst.eotest.image_utils as imutils
from lsst.eotest.sensor import MaskedCCD, add_mask_files, BrightPixels, \
    AmplifierGeometry, generate_mask
import lsst.eotest.sensor.sim_tools as sim_tools
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath

class BiasFunc(object):
    def __init__(self, slope, intercept):
        self.pars = slope, intercept
    def __call__(self, x):
        return x*self.pars[0] + self.pars[1]

class MaskedCCDTestCase(unittest.TestCase):
    gain = 1
    exptime = 1
    signal = 1
    xmin, xmax = 200, 250
    ymin, ymax = 1000, 1050
    mask_image = 'mask_image.fits'
    mpd = dict(afwImage.Mask().getMaskPlaneDict().items())
    @classmethod
    def setUpClass(cls):
        ccd = sim_tools.CCD(exptime=cls.exptime, gain=cls.gain)
        for amp in ccd.segments:
            imarr = ccd.segments[amp].image.getArray()
            imarr[cls.ymin:cls.ymax, cls.xmin:cls.xmax] += cls.signal
        ccd.writeto(cls.mask_image)
        cls.mask_files = []
        for mask_plane, bit in cls.mpd.items():
            mask_file = 'mask_file_%s.fits' % mask_plane
            cls.mask_files.append(mask_file)
            masked_ccd = MaskedCCD(cls.mask_image)
            pixels, columns = {}, {}
            for amp in masked_ccd:
                bp = BrightPixels(masked_ccd, amp, cls.exptime,
                                  cls.gain, mask_plane=mask_plane,
                                  ethresh=cls.signal/2.)
                pixels[amp], columns[amp] = bp.find()
            generate_mask(cls.mask_image, mask_file, mask_plane,
                          pixels=pixels, columns=columns)
        cls.summed_mask_file = 'summed_mask_file.fits'
        add_mask_files(cls.mask_files, cls.summed_mask_file)
    @classmethod
    def tearDownClass(cls):
        os.remove(cls.mask_image)
        for mask_file in cls.mask_files:
            os.remove(mask_file)
        os.remove(cls.summed_mask_file)
#    @unittest.skip('skip test_add_masks')
    def test_add_masks(self):
        ccd = MaskedCCD(self.mask_image)
        ccd.add_masks(self.summed_mask_file)
        total_signal = (self.signal*(self.ymax - self.ymin)
                        *(self.xmax - self.xmin))
        ny, nx = ccd[ccd.keys()[0]].getImage().getArray().shape
        for amp in ccd:
            self.assertEqual(sum(ccd[amp].getImage().getArray().flat),
                             total_signal)
        stat_ctrl = ccd.setMask('BAD')
        for amp in ccd:
            stats = afwMath.makeStatistics(ccd[amp], afwMath.MEAN, stat_ctrl)
            self.assertEqual(0, stats.getValue(afwMath.MEAN))
        stat_ctrl = ccd.setMask(clear=True)
        for amp in ccd:
            stats = afwMath.makeStatistics(ccd[amp], afwMath.MEAN, stat_ctrl)
            self.assertAlmostEqual(float(total_signal)/(nx*ny),
                                   stats.getValue(afwMath.MEAN), places=10)
#    @unittest.skip('skip test_setMask')
    def test_setMask(self):
        ccd = MaskedCCD(self.mask_image)
        for mp, bit in self.mpd.items():
            sctrl = ccd.setMask(mask_name=mp, clear=True)
            self.assertEqual(sctrl.getAndMask(), 2**bit)

class MaskedCCD_biasHandlingTestCase(unittest.TestCase):
    bias_slope = 1e-3
    bias_intercept = 0.5
    exptime = 1
    gain = 1
    image_file = 'test_image.fits'
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
        ccd = sim_tools.CCD(exptime=cls.exptime, gain=cls.gain)
        for amp in ccd.segments:
            ccd.segments[amp].image += cls.bias_image
        ccd.writeto(cls.image_file)

    @classmethod
    def tearDownClass(cls):
        os.remove(cls.image_file)

#    @unittest.skip('skip test_bias_image')
    def test_bias_image(self):
        ccd = MaskedCCD(self.image_file)
        for amp in ccd:
            my_bias_image = ccd.bias_image(amp)
            fracdiff = ( (self.bias_image.getArray()-my_bias_image.getArray())
                         /self.bias_image.getArray() )
            self.assertTrue(max(np.abs(fracdiff.flat)) < 1e-6)

#    @unittest.skip('skip test_unbias_and_trim')
    def test_unbias_and_trim(self):
        ccd = MaskedCCD(self.image_file)
        for amp in ccd:
            #
            # Test of corresponding MaskedCCD method.
            #
            image = ccd.unbiased_and_trimmed_image(amp)
            imarr = image.getImage().getArray()
            self.assertTrue(max(np.abs(imarr.flat)) < 1e-6)

    def test_bias_frame_subtraction(self):
        bias_file = 'bias_image.fits'
        bias_func_x = lambda x: 100*(np.sin(x/509.*2.*np.pi) + 1)
        bias_func_y = BiasFunc(self.bias_slope, self.bias_intercept)
        bias_ccd = sim_tools.CCD(exptime=self.exptime, gain=self.gain)
        for amp in bias_ccd.segments:
            imarr = bias_ccd.segments[amp].image.getArray()
            yvals = np.arange(0, imarr.shape[0], dtype=np.float)
            for x in range(imarr.shape[1]):
                imarr[:, x] += bias_func_x(x) + bias_func_y(yvals)
        bias_ccd.writeto(bias_file)

        image_file = 'image_file.fits'
        signal_level = 1000.
        image_ccd = sim_tools.CCD(exptime=self.exptime, gain=self.gain)
        for amp in image_ccd.segments:
            imarr = image_ccd.segments[amp].image.getArray()
            yvals = np.arange(0, imarr.shape[0], dtype=np.float)
            for x in range(imarr.shape[1]):
                imarr[:, x] +=\
                    bias_func_x(x) + bias_func_y(yvals)
            imaging_section = image_ccd.segments[amp].image.Factory(image_ccd.segments[amp].image, self.amp_geom.imaging)
            imaging_section += signal_level
        image_ccd.writeto(image_file)

        ccd = MaskedCCD(image_file, bias_frame=bias_file)
        for amp in ccd:
            image = ccd.unbiased_and_trimmed_image(amp)
            self.assertAlmostEqual(max(np.abs(image.getImage().getArray().ravel()))/signal_level, 1, 6)

        ccd = MaskedCCD(image_file)
        for amp in ccd:
            image = ccd.unbiased_and_trimmed_image(amp)
            self.assertTrue(max(np.abs(image.getImage().getArray().ravel()))/signal_level > 1.01)

        os.remove(bias_file)
        os.remove(image_file)

if __name__ == '__main__':
    unittest.main()
