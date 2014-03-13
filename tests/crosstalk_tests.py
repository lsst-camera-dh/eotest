"""
@brief Unit tests for crosstalk module.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import unittest
import lsst.eotest.image_utils as imutils
try:
    from lsst.eotest.sensor import MaskedCCD
    import lsst.eotest.sensor.sim_tools as sim_tools
    import lsst.eotest.sensor.crosstalk as crosstalk
except ImportError:
    # This is to allow this unit test to run on the inadequately
    # configured lsst-build01 on which Jenkins at SLAC runs.
    print "Error importing lsst.eotest.sensor"
    import sys
    sys.path.insert(0, os.path.join(os.environ['TEST_SCRIPTS_DIR'],
                                    'python', 'lsst', 'eotest', 'sensor'))
    from MaskedCCD import MaskedCCD
    import sim_tools
    import crosstalk

class CrosstalkTestCase(unittest.TestCase):
    """Test case for crosstalk code."""
    def setUp(self):
        self.xtalk_file = 'xtalk_test.fits'
        self.aggressor = 6
        dn = 2000
        x, y, radius = 250, 250, 20
        self.xtalk_frac = sim_tools.xtalk_pattern(self.aggressor)
        ccd = generate_crosstalk_frame(self.aggressor, dn, x, y, radius,
                                       xtalk_frac=self.xtalk_frac)
        ccd.writeto(self.xtalk_file)
    def tearDown(self):
        os.remove(self.xtalk_file)
    def test_detector_crosstalk(self):
        ccd = MaskedCCD(self.xtalk_file)
        ratios = crosstalk.detector_crosstalk(ccd, self.aggressor)
        for amp in ratios:
            if amp != self.aggressor:
                self.assertTrue(abs(ratios[amp][0] - self.xtalk_frac[amp])
                                < ratios[amp][1])

class CrosstalkMatrixTestCase(unittest.TestCase):
    def setUp(self):
        self.matrix_output = 'xtalk_output.txt'
        self.xtalk_files = []
        for agg in imutils.allAmps:
            self.xtalk_files.append('xtalk_test_%02i.fits' % agg)
            xtalk_frac = sim_tools.xtalk_pattern(agg)
            ccd = generate_crosstalk_frame(agg, 2000, 250, 250, 20,
                                           xtalk_frac=xtalk_frac)
            ccd.writeto(self.xtalk_files[-1])
    def tearDown(self):
        for item in self.xtalk_files:
            os.remove(item)
        os.remove(self.matrix_output)
        os.remove(self.matrix_output.replace('.txt', '.fits'))
    def test_CrosstalkMatrix(self):
        det_xtalk = crosstalk.make_crosstalk_matrix(self.xtalk_files,
                                                    verbose=False)
        det_xtalk.write(self.matrix_output)
        det_xtalk.write_fits(self.matrix_output.replace('.txt', '.fits'))

        det_xtalk2 = crosstalk.CrosstalkMatrix(self.matrix_output)
        diff = det_xtalk - det_xtalk2
        self.assertTrue(max([abs(x) for x in diff.matrix.flat]) < 1e-4)

def generate_crosstalk_frame(aggressor, dn, x, y, radius,
                             xtalk_frac=None, nom_frac=0.1):
    if xtalk_frac is None:
        xtalk_frac = dict([(amp, nom_frac) for amp in imutils.allAmps])
    ccd = sim_tools.CCD()
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
