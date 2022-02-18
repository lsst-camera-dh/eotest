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
    def run(self, sensor_id, detrespfile, gains,
            fit_range=(0., 9e4),
            nprofile_bins=100,
            unit_point=None,
            outputfile=None,
            plotfile=None,
            spline_ext_method=0,
            spline_s_factor=None):
        """Run the analysis

        Parameters
        ----------
        sensor_id : `str`
            The name of the CCD we are analyzing
        detrespfile : `str`
            The detector response file, produced by `FlatPairTask`
        gains : `dict`
            Dictionary of gains keyed by amp number.
        fit_range : `tuple`
            The range (in e-) over which to fit the splines
        nprofile_bins : `int`
            The number of bin to use when making the profile histogram
        null_point : `float`
            The ADU value at which the nonlinearity correction should vanish
        outputfile : `str` or `None`
            The name of the file to write the nonlinearity correction to
        plotfile : `str` or `None`
            The name of the file to write the diagnostic plots to
        spline_ext_method : `int` or `None`
            The method to use to extrapolate the spline
        spline_s_factor : `float` or `None`

        This make a profile plot of the fractional residual, (mean - q*slope) / q*slope,
        of the photodiode flux version amplifier mean for each amplifier, taking the
        values from the detrespfile DetectorResponse file.

        It then store the profile values in the outputfile FITS file.

        These value can then be used to construct a spline or implement whatver other correction.
        By default the `NonlinearityCorrection` class will construct a 3rd order spline for each amp
        """
        import sys
        self.sensor_id = sensor_id
        self.detrespfile = detrespfile
        self.fit_range = fit_range
        self.nprofile_bins = nprofile_bins
        self.null_point = null_point
        self.detresp = DetectorResponse(detrespfile)

        kw_ctor = {}
        if spline_ext_method is not None:
            kw_ctor['ext'] = spline_ext_method
        if spline_s_factor is not None:
            kw_ctor['s'] = spline_s_factor

        self.nlc = NonlinearityCorrection.create_from_det_response(self.detresp,
                                                                   gains,
                                                                   fit_range=self.fit_range,
                                                                   nprofile_bins=self.nprofile_bins,
                                                                   null_point=self.null_point,
                                                                   **kw_ctor)

        if outputfile is not None:
            fulloutput = os.path.join(self.config.output_dir, outputfile)
            self.nlc.write_to_fits(fulloutput)

        if plotfile is not None:
            fullplotpath = os.path.join(self.config.output_dir, plotfile)
            self.nlc.save_plots(fullplotpath)

