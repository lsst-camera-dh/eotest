"""
@brief Unit tests for ccd250_mask module.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import unittest
import numpy as np
import pyfits
import lsst.test_scripts.image_utils as imutils
from lsst.test_scripts.sensor.ccd250_mask import ccd250_mask

class _FitsFile(dict):
    def __init__(self, infile):
        dict.__init__(self)
        foo = pyfits.open(infile)
        for amp in imutils.allAmps:
            self[amp] = foo[amp].data

class ccd250_mask_TestCase(unittest.TestCase):
    """Test case for ccd250_mask function."""
    def setUp(self):
        self.mask_file = 'ccd250_defects_mask.fits'
        self.image_file = 'temp_mask_image.fits'
        self.outer_edge_width = 10
        self.bloom_stop_width = 5
        self.signal = 10
    def tearDown(self):
        os.remove(self.mask_file)
        os.remove(self.image_file)
    def test_ccd250_mask(self):
        ccd250_mask(self.mask_file, tmp_mask_image=self.image_file,
                    outer_edge_width=self.outer_edge_width,
                    bloom_stop_width=self.bloom_stop_width,
                    signal=self.signal, cleanup=False)
        image = _FitsFile(self.image_file)
        mask = _FitsFile(self.mask_file)
        for amp in imutils.allAmps:
            #
            # Unmasked region.
            #
            indx = np.where(image[amp] == 0)
            #
            # Verify expected sensor perimeter mask along vertical sides.
            #
            self.assertEqual(min(indx[0]),
                             imutils.imaging.getMinY() + self.outer_edge_width)
            #
            # Verify that mask has zero bits set in unmasked region.
            #
            self.assertEqual(min(mask[amp][indx].flat), 0)
            self.assertEqual(max(mask[amp][indx].flat), 0)
            #
            # Check that mask area is subset of mask image area.
            #
            indx = np.where(mask[amp] != 0)
            self.assertTrue(min(image[amp][indx].flat) >= self.signal)

if __name__ == '__main__':
    unittest.main()
