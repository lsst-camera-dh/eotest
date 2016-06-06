"""
@brief Unit tests for rolloff_mask module.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import unittest
import numpy as np
import astropy.io.fits as fits
import lsst.eotest.image_utils as imutils
from lsst.eotest.sensor.rolloff_mask import rolloff_mask
from lsst.eotest.sensor import AmplifierGeometry, makeAmplifierGeometry, amp_loc
import lsst.eotest.sensor.sim_tools as sim_tools

class _FitsFile(dict):
    def __init__(self, infile):
        dict.__init__(self)
        foo = fits.open(infile)
        amps = imutils.allAmps(infile)
        for amp in amps:
            self[amp] = foo[amp].data

class rolloff_mask_TestCase(unittest.TestCase):
    """Test case for rolloff_mask function."""
    def setUp(self):
        self.input_file = 'input_image.fits'
        self.mask_file = 'rolloff_defects_mask.fits'
        self.image_file = 'my_temp_mask_image.fits'
        self.outer_edge_width = 10
        self.bloom_stop_width = 5
        self.signal = 10
        ccd = sim_tools.CCD(geometry=AmplifierGeometry(amp_loc=amp_loc['E2V']))
        ccd.writeto(self.input_file, bitpix=16)
    def tearDown(self):
        os.remove(self.mask_file)
        os.remove(self.image_file)
        os.remove(self.input_file)
    def test_rolloff_mask(self):
        amp_geom = makeAmplifierGeometry(self.input_file)
        rolloff_mask(self.input_file, self.mask_file,
                     tmp_mask_image=self.image_file,
                     outer_edge_width=self.outer_edge_width,
                     bloom_stop_width=self.bloom_stop_width,
                     signal=self.signal, cleanup=False)
        image = _FitsFile(self.image_file)
        mask = _FitsFile(self.mask_file)
        for amp in imutils.allAmps(self.input_file):
            #
            # Unmasked region.
            #
            indx = np.where(image[amp] == 0)
            #
            # Verify expected sensor perimeter mask along vertical sides.
            #
            self.assertEqual(min(indx[0]),
                             amp_geom.imaging.getMinY() + self.outer_edge_width)
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
