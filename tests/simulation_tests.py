"""
@brief Tests of simulation.sim_tools.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import unittest
import numpy as np
from simulation.sim_tools import *

class SegmentExposureTestCase(unittest.TestCase):
    def setUp(self):
        self.seg = SegmentExposure(exptime=100, gain=1, ccdtemp=-100)
        self.intensity = 100
        self.full_well = 1.5e5
    def tearDown(self):
        pass
    def test_expose_flat(self):
        times = np.arange(0, 1000, self.seg.exptime)
        for i, time in enumerate(times):
            image = imutils.unbias_and_trim(self.seg.image)
            image_mean = imutils.mean(image)
            illum = i*self.intensity*self.seg.exptime/self.seg.gain
            if i != 0:
                self.assertTrue((illum - image_mean)/illum < 3e-4)
            self.seg.expose_flat(intensity=self.intensity)
    def test_full_well(self):
        self.seg.full_well = self.full_well
        times = np.arange(0, 2000, self.seg.exptime)
        for i, time in enumerate(times):
            image = imutils.unbias_and_trim(self.seg.image)
            Ne_mean = imutils.mean(image)*self.seg.gain
            self.assertTrue(Ne_mean <= self.full_well)
            self.seg.expose_flat(intensity=self.intensity)
        
if __name__ == '__main__':
    unittest.main()
