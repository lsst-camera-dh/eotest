"""
@brief Unit tests for EOTestResults.py module.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
from builtins import range
import os
import unittest
import numpy.random as ra
from lsst.eotest.sensor import EOTestResults


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

    def test_append_column(self):
        results = EOTestResults(self.results_file)
        results.append_column('NEW_COLUMN', int)
        values = results['NEW_COLUMN']
        for i in range(16):
            self.assertEqual(values[i], 0)

    def test_add_seg_result_new_column(self):
        results = EOTestResults(self.results_file)
        for amp in range(1, 17):
            results.add_seg_result(amp, 'NEW_COLUMN', 2.1)
        values = results['NEW_COLUMN']
        for i in range(16):
            self.assertAlmostEqual(values[i], 2.1, places=6)


if __name__ == '__main__':
    unittest.main()
