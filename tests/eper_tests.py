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
    Base class for eper TestCase classes.
    """
    imaging_value = 100
    overscan_value = 1
    overscans = 2
    verbose = False


class SerialEperTestCase(EperTestCase):
    """
    TestCase class for EPERTask implementation in the serial direction.
    """

    def setUp(self):
        """
        Create a CCD frame FITS file with known imaging and overscan
        values.
        """
        self.fits_file = 'eper_serial_test_frame.fits'
        ccd = sim_tools.CCD()
        for amp in ccd.segments:
            segment = ccd.segments[amp]
            segment.imarr[:, -1] = self.imaging_value
            overscan_region = segment.geometry.serial_overscan
            overscan = segment.image.Factory(segment.image,
                                             overscan_region).getArray()
            overscan[:, 0:self.overscans] = self.overscan_value
        ccd.writeto(self.fits_file)

    def tearDown(self):
        "Delete test FITS file."
        os.remove(self.fits_file)

    def test_subimage_serial(self):
        "Test the SubImage class for the serial direction."
        ccd = sensorTest.MaskedCCD(self.fits_file)
        task = sensorTest.EPERTask()
        task.config.direction = 's'
        task.config.verbose = self.verbose
        for amp in ccd:
            subimage = SubImage(ccd, amp, self.overscans, task)
            lastcol = ccd.amp_geom.imaging.getMaxX()
            nrows = ccd.amp_geom.imaging.getHeight()
            last_imaging_col = subimage(lastcol).getImage().getArray().flatten()
            self.assertEqual(nrows, len(last_imaging_col))
            self.assertEqual(nrows*self.imaging_value, sum(last_imaging_col))
            for icol in range(1, self.overscans+1):
                overscan_col = subimage(lastcol + icol)
                self.assertEqual(nrows*self.overscan_value,
                                 sum(overscan_col.getImage().getArray().flatten()))

    def test_EPERTask_run_serial(self):
        "Test the EPERTask.run method for the serial direction."
        ccd = sensorTest.MaskedCCD(self.fits_file)
        task = sensorTest.EPERTask()
        task.config.direction = 's'
        task.config.cti = True
        task.config.verbose = self.verbose
        cti, bias_ests = task.run(self.fits_file, range(1, 17), self.overscans)
        ncols = (ccd.amp_geom.prescan.getWidth() +
                 ccd.amp_geom.imaging.getWidth())
        cti_expected = (float(self.overscans*self.overscan_value)
                        / float(self.imaging_value)/float(ncols))
        for amp in ccd:
            self.assertEqual(cti[amp].value, cti_expected)


class ParallelEperTestCase(EperTestCase):
    """
    TestCase class for EPERTask implementation in the parallel direction.
    """

    def setUp(self):
        """
        Create a CCD frame FITS file with known imaging and overscan
        values.
        """
        self.fits_file = 'eper_parallel_test_frame.fits'
        ccd = sim_tools.CCD()
        for amp in ccd.segments:
            segment = ccd.segments[amp]
            segment.imarr[-1, :] = self.imaging_value
            overscan_region = segment.geometry.parallel_overscan
            overscan = segment.image.Factory(segment.image,
                                             overscan_region).getArray()
            overscan[0:self.overscans, :] = self.overscan_value
        ccd.writeto(self.fits_file)

    def tearDown(self):
        "Delete test FITS file."
        os.remove(self.fits_file)

    def test_subimage_parallel(self):
        "Test the SubImage class for the parallel direction."
        ccd = sensorTest.MaskedCCD(self.fits_file)
        task = sensorTest.EPERTask()
        task.config.direction = 'p'
        task.config.verbose = self.verbose
        for amp in ccd:
            subimage = SubImage(ccd, amp, self.overscans, task)
            lastrow = ccd.amp_geom.imaging.getMaxY()
            ncols = ccd.amp_geom.imaging.getWidth()
            last_imaging_row = subimage(lastrow).getImage().getArray().flatten()
            self.assertEqual(ncols, len(last_imaging_row))
            self.assertEqual(ncols*self.imaging_value, sum(last_imaging_row))
            for irow in range(1, self.overscans+1):
                overscan_row = subimage(lastrow + irow)
                self.assertEqual(ncols*self.overscan_value,
                                 sum(overscan_row.getImage().getArray().flatten()))

    def test_EPERTask_run_parallel(self):
        "Test the EPERTask.run method for the parallel direction."
        ccd = sensorTest.MaskedCCD(self.fits_file)
        task = sensorTest.EPERTask()
        task.config.direction = 'p'
        task.config.cti = True
        task.config.verbose = self.verbose
        cti, bias_ests = task.run(self.fits_file, range(1, 17), self.overscans)
        nrows = ccd.amp_geom.imaging.getHeight()
        cti_expected = (float(self.overscans*self.overscan_value)
                        / float(self.imaging_value)/float(nrows))
        for amp in ccd:
            self.assertEqual(cti[amp].value, cti_expected)


if __name__ == '__main__':
    unittest.main()
