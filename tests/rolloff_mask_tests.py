"""
@brief Unit tests for rolloff_mask module.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import copy
import unittest
import numpy as np
import astropy.io.fits as fits
import lsst.eotest.image_utils as imutils
from lsst.eotest.sensor.rolloff_mask import rolloff_mask
from lsst.eotest.sensor import AmplifierGeometry, makeAmplifierGeometry, amp_loc
import lsst.eotest.sensor.sim_tools as sim_tools


class _FitsFile(dict):
    def __init__(self, infile):
        amp_geom = makeAmplifierGeometry(infile)
        xmin, xmax = amp_geom.imaging.getMinX(), amp_geom.imaging.getMaxX()
        super(_FitsFile, self).__init__()
        with fits.open(infile) as foo:
            amps = imutils.allAmps(infile)
            for amp in amps:
                self[amp] = copy.deepcopy(foo[amp].data[:, xmin:xmax])


class rolloff_mask_TestCase(unittest.TestCase):
    """Test case for rolloff_mask function."""

    def setUp(self):
        self.input_file = 'input_image.fits'
        self.mask_file = 'rolloff_defects_mask.fits'
        self.image_file = 'my_temp_mask_image.fits'
        self.outer_edge_width = 10
        self.bloom_stop_width = 5
        self.signal = 10
        self.amp_geoms \
            = (AmplifierGeometry(amp_loc=amp_loc['E2V']),
               AmplifierGeometry(prescan=3, nx=509, ny=2000,
                                 amp_loc=amp_loc['ITL'], vendor='ITL'))

    def tearDown(self):
        os.remove(self.mask_file)
        os.remove(self.image_file)
        os.remove(self.input_file)

    def test_rolloff_mask(self):
        for amp_geom in self.amp_geoms:
            ccd = sim_tools.CCD(geometry=amp_geom)
            ccd.writeto(self.input_file, bitpix=16)
            rolloff_mask(self.input_file, self.mask_file,
                         tmp_mask_image=self.image_file,
                         outer_edge_width=self.outer_edge_width,
                         bloom_stop_width=self.bloom_stop_width,
                         signal=self.signal, cleanup=False)
            image = _FitsFile(self.image_file)
            mask = _FitsFile(self.mask_file)
            ymin = amp_geom.imaging.getMinY()
            ymax = amp_geom.imaging.getMaxY()
            for amp in imutils.allAmps(self.input_file):
                #
                # Unmasked region.
                #
                indx = np.where(image[amp] == 0)
                #
                # Verify expected sensor perimeter mask along vertical sides.
                #
                self.assertEqual(min(indx[0]), ymin + self.outer_edge_width)
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
            # Check masked regions along sensor edges for amps 1, 8, 9, 16:
            if amp_geom.vendor == 'E2V':
                for amp in (1, 9):
                    self.assertGreaterEqual(image[amp][ymax, -1], self.signal)
                for amp in (8, 16):
                    self.assertGreaterEqual(image[amp][ymax, 0], self.signal)
            if amp_geom.vendor == 'ITL':
                for amp in (1, 16):
                    self.assertGreaterEqual(image[amp][ymax, -1], self.signal)
                for amp in (8, 9):
                    self.assertGreaterEqual(image[amp][ymax, 0], self.signal)

if __name__ == '__main__':
    unittest.main()
