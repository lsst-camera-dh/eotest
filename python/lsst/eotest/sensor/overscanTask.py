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
from .overscan import OverscanResults

class OverscanConfig(pexConfig.Config):
    """Configuration for overscan analysis task"""
    output_dir = pexConfig.Field("Output directory", str, default='.')
    output_file = pexConfig.Field("Output filename", str, default=None)
    verbose = pexConfig.Field("Turn verbosity on", bool, default=True)

class OverscanTask(pipeBase.Task):
    """Task to calculate mean row results from flat pair dataset."""
    ConfigClass = OverscanConfig
    _DefaultName = "OverscanTask"

    def run(self, sensor_id, infiles, gains, bias_frame=None,
            linearity_correction=None):

        all_amps = imutils.allAmps(infiles[0])
        ## Calculate mean row for each flat file
        overscan_results = OverscanResults(all_amps)
        for i, infile in enumerate(infiles):
            if self.config.verbose:
                self.log.info("Processing {0}".format(infile))
            ccd = MaskedCCD(infile, bias_frame=bias_frame,
                            linearity_correction=linearity_correction)
            overscan_results.process_image(ccd, gains)
                
        output_dir = self.config.output_dir
        if self.config.output_file is None:
            output_file = os.path.join(output_dir, 
                                       '{0}_overscan_results.fits'.format(sensor_id))
        else:
            output_file = os.path.join(output_dir, self.config.output_file)
        if self.config.verbose:
            self.log.info("writing to {0}".format(output_file))
        overscan_results.write_results(outfile=output_file)
