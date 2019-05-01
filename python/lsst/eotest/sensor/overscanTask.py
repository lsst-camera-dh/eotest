"""
@brief Additional overscan analysis on flat pair eotest
"""
from __future__ import print_function
from __future__ import absolute_import
import os
import numpy as np
import astropy.io.fits as fits

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.eotest.image_utils as imutils
from lsst.eotest.fitsTools import fitsWriteto
from lsst.eotest.sensor import MaskedCCD, parse_geom_kwd
from .overscan_fit import OverscanFit

class OverscanConfig(pexConfig.Config):
    """Configuration for overscan analysis task"""
    output_dir = pexConfig.Field("Output directory", str, default='.')
    output_file = pexConfig.Field("Output filename", str, default=None)
    smoothing = pexConfig.Field("Smoothing for spline overscan correction",
                                int, default=11000)
    minflux = pexConfig.Field("Minimum flux for overscan fitting.", float,
                              default=10000.0)
    maxflux = pexConfig.Field("Maximum flux for overscan fitting.", float,
                              default=140000.0)
    num_oscan_pixels = pexConfig.Field("Number of overscan pixels used for model fit.",
                                       int, default=10)
    verbose = pexConfig.Field("Turn verbosity on", bool, default=True)

class OverscanTask(pipeBase.Task):
    """Task to calculate mean row results from flat pair dataset."""
    ConfigClass = OverscanConfig
    _DefaultName = "OverscanTask"

    def run(self, sensor_id, infiles, gains, bias_frame=None):

        ## Calculate mean row for each flat file
        minflux = self.config.minflux
        maxflux = self.config.maxflux
        num_oscan_pixels = self.config.num_oscan_pixels

        fitter = OverscanFit(num_oscan_pixels=num_oscan_pixels, minflux=minflux, maxflux=maxflux)
        for i, infile in enumerate(infiles):
            if self.config.verbose:
                self.log.info("Processing {0}".format(infile))
            ccd = MaskedCCD(infile, bias_frame=bias_frame)
            fitter.process_image(ccd, gains)
                
        output_dir = self.config.output_dir
        if self.config.output_file is None:
            output_file = os.path.join(output_dir, 
                                       '{0}_overscan_results.fits'.format(sensor_id))
        else:
            output_file = os.path.join(output_dir, self.config.output_file)
        if self.config.verbose:
            self.log.info("writing to {0}".format(output_file))
        fitter.write_results(outfile=output_file)

        return output_file
