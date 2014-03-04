"""
@brief Unit tests for EOTestResults.py module.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import unittest
import numpy.random as ra

try:
    from lsst.eotest.sensor import EOTestResults
except ImportError:
    # This is to allow this unit test to run on the inadequately
    # configured lsst-build01 on which Jenkins at SLAC runs.
    print "Error importing lsst.eotest.sensor"
    import sys
    sys.path.insert(0, os.path.join(os.environ['TEST_SCRIPTS_DIR'],
                                    'python', 'lsst', 'eotest', 'sensor'))
    from EOTestResults import EOTestResults

class EOTestResultsTestCase(unittest.TestCase):
    def setUp(self):
        self.results_file = 'test_results.fits'
        results = EOTestResults(self.results_file)
        self.gains = [ra.normal(0, 1) for x in range(16)]
        for amp in range(1, 17):
            results.add_seg_result(amp, 'GAIN', self.gains[amp-1])
        results.write(clobber=True)
    def tearDown(self):
        try:
            os.remove(self.results_file)
        except OSError:
            pass
    def test_retrieve_values(self):
        results = EOTestResults(self.results_file)
        gains = results['GAIN']
        for i in range(16):
            self.assertAlmostEqual(gains[i], self.gains[i], places=6)

if __name__ == '__main__':
    unittest.main()


