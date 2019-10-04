"""
@brief Pixel response non-linearity.

@author E. Charles <echarles@slac.stanford.edu>
"""
from __future__ import absolute_import

import os

import lsst.pipe.base as pipeBase
import lsst.pex.config as pexConfig

from lsst.eotest.sensor.DetectorResponse import DetectorResponse
from lsst.eotest.sensor.NonlinearityCorrection import NonlinearityCorrection

class NonlinearityConfig(pexConfig.Config):
    """Configuration for pixel response non-uniformity task"""
    output_dir = pexConfig.Field("Output directory", str, default=".")
    verbose = pexConfig.Field("Turn verbosity on", bool, default=True)

class NonlinearityTask(pipeBase.Task):
    """Task for computing pixel nonlinearity"""
    ConfigClass = NonlinearityConfig
    _DefaultName = "NonlinearityTask"

    @pipeBase.timeMethod
    def run(self, sensor_id, detrespfile,
            fit_range=(0., 9e4),
            nprofile_bins=10,
            outputfile=None,
            plotfile=None):
        """Run the analysis

        Parameters
        ----------
        sensor_id : `str`
            The name of the CCD we are analyzing
        detrespfile : `str`
            The detector response file, produced by `FlatPairTask`
        fit_range : `tuple`
            The range (in e-) over which to fit the splines
        nprofile_bins : `int`
            The number of bin to use when making the profile histogram
        outputfile : `str` or `None`
            The name of the file to write the nonlinearity correction to
        plotfile : `str` or `None`
            The name of the file to write the diagnostic plots to

        This make a profile plot of the fractional residual, (mean - q*slope) / q*slope,
        of the photodiode flux version amplifier mean for each amplifier, taking the
        values from the detrespfile DetectorResponse file.

        It then store the profile values in the outputfile FITS file.

        These value can then be used to construct a spline or implement whatver other correction.
        By default the `NonlinearityCorrection` class will construct a 3rd order spline for each amp
        """

        self.sensor_id = sensor_id
        self.detrespfile = detrespfile
        self.fit_range = fit_range
        self.nprofile_bins = nprofile_bins
        self.detresp = DetectorResponse(detrespfile)

        self.nlc = NonlinearityCorrection.create_from_det_response(self.detresp,
                                                                   self.fit_range,
                                                                   self.nprofile_bins)

        if outputfile is not None:
            fulloutput = os.path.join(self.config.output_dir, outputfile)
            self.nlc.write_to_fits(fulloutput)

        if plotfile is not None:
            fullplotpath = os.path.join(self.config.output_dir, plotfile)
            self.nlc.save_plots(fullplotpath)



if __name__ == '__main__':

    import lsst.afw.math as afwMath
    from lsst.eotest.sensor import MaskedCCD

    sensorid = 'RTM-022'
    detrespfile1 = '/gpfs/slac/lsst/fs3/g/data/jobHarness/jh_archive/LCA-11021_RTM/LCA-11021_RTM-022/11671/flat_pairs_raft_analysis/v0/90861/S00/ITL-3800C-080_det_response.fits'

    bias_frame = '/gpfs/slac/lsst/fs3/g/data/jobHarness/jh_archive/LCA-11021_RTM/LCA-11021_RTM-022/11671/fe55_raft_analysis/v0/90848/S00/ITL-3800C-080_mean_bias_5.fits'
    flat_frame = '/gpfs/slac/lsst/fs3/g/data/jobHarness/jh_archive/LCA-11021_RTM/LCA-11021_RTM-022/11671/flat_pair_raft_acq/v0/90850/S00/ITL-3800C-080_flat_0000.09s_flat1_11671_20190824200500.fits'

    test_out = 'nonlin.fits'
    test_plot = 'nonlin_plot.png'

    task = NonlinearityTask()
    task.run(sensorid, detrespfile1, outputfile=test_out, plotfile=test_plot)

    nlc = NonlinearityCorrection.create_from_fits_file(test_out)

    ccd_1 = MaskedCCD(flat_frame, bias_frame=bias_frame)
    ccd_2 = MaskedCCD(flat_frame, bias_frame=bias_frame, linearity_correction=nlc)

    img_1 = ccd_1.unbiased_and_trimmed_image(1)
    img_2 = ccd_2.unbiased_and_trimmed_image(1)

    mean_1 = afwMath.makeStatistics(img_1, afwMath.MEAN, ccd_1.stat_ctrl).getValue()
    mean_2 = afwMath.makeStatistics(img_2, afwMath.MEAN, ccd_2.stat_ctrl).getValue()

    print(mean_1, mean_2)
