"""
Unit tests for extended pixel edge response code.
"""
import os
import unittest
import lsst.eotest.sensor as sensorTest

import lsst.afw.math as afwMath
from lsst.eotest.sensor import NonlinearityTask, MaskedCCD, NonlinearityCorrection 

_root_dir = '/gpfs/slac/lsst/fs3/g/data/jobHarness/jh_archive/LCA-11021_RTM/LCA-11021_RTM-022/11671/'
_detrespfile = os.path.join(_root_dir, 'flat_pairs_raft_analysis/v0/90861/S00/ITL-3800C-080_det_response.fits')
_bias_frame = os.path.join(_root_dir, 'fe55_raft_analysis/v0/90848/S00/ITL-3800C-080_mean_bias_5.fits')
_flat_frame = os.path.join(_root_dir, 'flat_pair_raft_acq/v0/90850/S00/ITL-3800C-080_flat_0000.09s_flat1_11671_20190824200500.fits')


class NonlinearityCorrectionTestCase(unittest.TestCase):
    """
    TestCase class for NonlinearityCorrection
    """

    def setUp(self):
        """
        Just set the name of the output file we will be making
        """
        self.fits_file = 'nonlin.fits'

    def tearDown(self):
        "Delete test FITS file."
        os.remove(self.fits_file)


    @unittest.skipUnless(os.path.isdir(_root_dir), 'SLAC data files not avaiable.')
    def test_01_nonlinearity_task(self):
        """Run the NonlinearityTask"""
        task = NonlinearityTask()
        sensorid = 'RTM-022'
        task.run(sensorid, _detrespfile, outputfile=self.fits_file, plotfile=None)

        nlc = NonlinearityCorrection.create_from_fits_file(self.fits_file)
        
        ccd_1 = MaskedCCD(_flat_frame, bias_frame=_bias_frame)
        ccd_2 = MaskedCCD(_flat_frame, bias_frame=_bias_frame, linearity_correction=nlc)

        img_1 = ccd_1.unbiased_and_trimmed_image(1)
        img_2 = ccd_2.unbiased_and_trimmed_image(1)
        
        mean_1 = afwMath.makeStatistics(img_1, afwMath.MEAN, ccd_1.stat_ctrl).getValue()
        mean_2 = afwMath.makeStatistics(img_2, afwMath.MEAN, ccd_2.stat_ctrl).getValue()

        self.assertAlmostEqual(mean_1, 40.724625956352185)
        self.assertAlmostEqual(mean_2, 40.724101151614555)


if __name__ == '__main__':
    unittest.main()
