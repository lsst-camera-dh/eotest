import os
import unittest
from simulation.sim_tools import CCD
from crosstalk import *

class CrosstalkTestCase(unittest.TestCase):
    """Test case for crosstalk code."""
    def setUp(self):
        self.xtalk_file = 'xtalk_test.fits'
        self.aggressor = 6
        dn = 2000
        x, y, radius = 250, 250, 20
        self.xtalk_frac = {}
        nside = len(imutils.allAmps)/2
        for victim in imutils.allAmps:
            if (victim != self.aggressor and
                (victim-1)/nside == (self.aggressor-1)/nside):
                dist = abs(victim - self.aggressor)
                self.xtalk_frac[victim] = 0.02/dist**2
            else:
                self.xtalk_frac[victim] = 0
        ccd = generate_crosstalk_frame(self.aggressor, dn, x, y, radius,
                                       xtalk_frac=self.xtalk_frac)
        ccd.writeto(self.xtalk_file)
    def tearDown(self):
        os.remove(self.xtalk_file)
    def test_detector_crosstalk(self):
        ccd = MaskedCCD(self.xtalk_file)
        ratios = detector_crosstalk(ccd, self.aggressor)
        for amp in ratios:
            if amp != self.aggressor:
                self.assertTrue(abs(ratios[amp][0] - self.xtalk_frac[amp])
                                < ratios[amp][1])

def generate_crosstalk_frame(aggressor, dn, x, y, radius,
                             xtalk_frac=None, nom_frac=0.1):
    if xtalk_frac is None:
        xtalk_frac = dict([(amp, nom_frac) for amp in imutils.allAmps])
    ccd = CCD()
    ccd.add_bias()
    ccd.add_dark_current()
    for amp in ccd.segments:
        if amp == aggressor:
            ccd.segments[amp].add_spot_image(dn, x, y, radius)
        else:
            ccd.segments[amp].add_spot_image(dn*xtalk_frac[amp], x, y, radius)
    return ccd

if __name__ == '__main__':
    unittest.main()
