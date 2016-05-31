"""
Unit tests for extended pixel edge response code.
"""
import os
import unittest
import lsst.eotest.sensor as sensorTest
from lsst.eotest.sensor.eperTask import SubImage
import lsst.eotest.sensor.sim_tools as sim_tools

class EperTestCase(unittest.TestCase):
    """
    TestCase class for EPERTask implementation.
    """
    def setUp(self):
        self.fits_file = 'eper_test_frame.fits'
        self.imaging_value = 100
        self.overscan_value = 1
        self.overscans = 2
        ccd = sim_tools.CCD()
        for amp in ccd.segments:
            segment = ccd.segments[amp]
            segment.imarr[:,-1] = self.imaging_value
            overscan = segment.image.Factory(segment.image,
                                             segment.geometry.serial_overscan).getArray()
            overscan[:,0:self.overscans] = self.overscan_value
        ccd.writeto(self.fits_file)

    def tearDown(self):
        os.remove(self.fits_file)
        pass

    def test_subimage(self):
        "Test the SubImage class."
        ccd = sensorTest.MaskedCCD(self.fits_file)
        task = sensorTest.EPERTask()
        task.config.direction = 's'
        amp = 1
        subimage = SubImage(ccd, amp, self.overscans, task)
        lastcol = ccd.amp_geom.imaging.getMaxX()
        nrows = ccd.amp_geom.imaging.getHeight()
        last_imaging_column = subimage(lastcol).getImage().getArray().flatten()
        self.assertEqual(nrows, len(last_imaging_column))
        self.assertEqual(nrows*self.imaging_value, sum(last_imaging_column))
        for icol in range(1, self.overscans+1):
            overscan_col = subimage(lastcol+icol)
            self.assertEqual(nrows*self.overscan_value,
                             sum(overscan_col.getImage().getArray().flatten()))

    def test_EPERTask_run(self):
        "Test the EPERTask.run method."
        ccd = sensorTest.MaskedCCD(self.fits_file)
        task = sensorTest.EPERTask()
        task.config.direction = 's'
        task.config.cti = True

        cti = task.run(self.fits_file, range(1, 17), self.overscans)

        ncols = (ccd.amp_geom.prescan.getWidth() +
                 ccd.amp_geom.imaging.getWidth())
        cti_expected = (float(self.overscans*self.overscan_value)
                        /float(self.imaging_value)/float(ncols))
        for amp in ccd:
            self.assertEqual(cti[amp].value, cti_expected)


if __name__ == '__main__':
    unittest.main()
