"""
Unit tests for generate_mask code.
"""
import os
import unittest
import lsst.eotest.image_utils as imutils
import lsst.eotest.sensor as sensorTest
import lsst.eotest.sensor.sim_tools as sim_tools


class GenerateMaskTestCase(unittest.TestCase):
    """
    TestCase class for mask generation.
    """

    def setUp(self):
        "Set the physical pixels and columns to mask."
        self.template_file = 'template_frame.fits'
        ccd = sim_tools.CCD()
        ccd.writeto(self.template_file)
        self.mask_file = 'test_mask_file.fits'
        # These are physical pixel coordinates in the imaging section,
        # i.e., prescan pixels are not included in the x-coordinate.
        self.pixels = dict([(1, [(200, 1000), (500, 300), (140, 1499)]),
                            (2, [(342, 6), (50, 2), (403, 11), (420, 12)])])
        self.columns = dict([(1, (120, 212, 320, 432))])

    def tearDown(self):
        "Clean up the temporary files."
        os.remove(self.template_file)
        os.remove(self.mask_file)

    def test_generate_mask(self):
        "Test the generated mask for expected masked and unmasked pixels."
        sensorTest.generate_mask(self.template_file,
                                 self.mask_file,
                                 mask_plane='TRAPS',
                                 pixels=self.pixels,
                                 columns=self.columns,
                                 temp_mask_image='my_temp_mask_file.fits')
        ccd = sensorTest.MaskedCCD(self.mask_file)
        for amp in self.pixels:
            image = imutils.trim(ccd[amp].getImage(), ccd.amp_geom.imaging)
            imarr = image.getArray()
            for ix, iy in self.pixels[amp]:
                self.assertNotEqual(0, imarr[iy][ix])
                for xoffset, yoffset in ((-1, -1), (-1, 0), (-1, 1), (0, -1),
                                         (0, 1), (1, -1), (1, 0), (1, 1)):
                    self.assertEqual(0, imarr[iy+yoffset][ix+xoffset])
        for amp in self.columns:
            image = imutils.trim(ccd[amp].getImage(), ccd.amp_geom.imaging)
            imarr = image.getArray()
            for ix in self.columns[amp]:
                self.assertNotEqual(0, imarr[0][ix])
                self.assertEqual(imarr[0][ix]*imarr.shape[0], sum(imarr[:, ix]))


if __name__ == '__main__':
    unittest.main()
