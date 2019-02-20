"""
@brief Source identification and characterization of
spot images.
"""
from __future__ import print_function
from __future__ import absolute_import
import os
import numpy as np
from astropy.io import fits

import lsst.eotest.image_utils as imutils
import lsst.afw.image as afwImage
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.tasks.characterizeImage import CharacterizeImageTask, CharacterizeImageConfig
import lsst.meas.extensions.shapeHSM

from .MaskedCCD import MaskedCCD
from .AmplifierGeometry import parse_geom_kwd

class SpotConfig(pexConfig.Config):
    """Configuration for Spot analysis task"""
    minpixels = pexConfig.Field("Minimum number of pixels above detection threshold", 
                                int, default=10)
    nsig = pexConfig.Field("Source footprint threshold in number of standard deviations.", 
                           float, default=10)
    temp_set_point = pexConfig.Field("Required temperature (C) set point",
                                     float, default=-95.)
    temp_set_point_tol = pexConfig.Field("Required temperature set point tolerance (degrees C)",
                                         float, default=1.)
    output_dir = pexConfig.Field("Output directory", str, default='.')
    output_file = pexConfig.Field("Output filename", str, default=None)
    verbose = pexConfig.Field("Turn verbosity on", bool, default=True)

class SpotTask(pipeBase.Task):
    """Task to estimate spot moments from spot projector data."""

    ConfigClass = SpotConfig
    _DefaultName = "SpotTask"

    @pipeBase.timeMethod
    def run(self, sensor_id, infile, gains, bias_frame=None):

        if self.config.verbose:
            self.log.info("Input file:")
            self.log.info("  {0}".format(infile))
        #
        # Process a CCD image mosaic
        #
        if self.config.verbose:
            self.log.info("processing {0}".format(infile))
        image = imutils.make_ccd_mosaic(infile, bias_frame=bias_frame, gains=gains)
        exposure = afwImage.ExposureF(image.getBBox())
        exposure.setImage(image)
        #
        # Set up characterize task configuration
        #
        nsig = self.config.nsig
        minpixels = self.config.minpixels
        charConfig = CharacterizeImageConfig()
        charConfig.doMeasurePsf = False
        charConfig.doApCorr = False
        charConfig.repair.doCosmicRay = False
        charConfig.detection.minPixels = minpixels
        charConfig.detection.background.binSize = 10
        charConfig.detection.thresholdType = "stdev"
        charConfig.detection.thresholdValue = nsig
        charConfig.measurement.plugins.names |= ["ext_shapeHSM_HsmSourceMoments"]
        charTask = CharacterizeImageTask(config=charConfig)
        result = charTask.characterize(exposure)
        src = result.sourceCat
        if self.config.verbose:
            self.log.info("Detected {0} objects".format(len(src)))
        #
        # Save catalog results to file
        #
        output_dir = self.config.output_dir
        if self.config.output_file is None:
            output_file = os.path.join(output_dir,
                                       '{0}_source_catalog.cat'.format(sensor_id))
        else:
            output_file = os.path.join(output_dir, self.config.output_file)
        if self.config.verbose:
            self.log.info("Writing spot results file to {0}".format(output_file))
        src.writeFits(output_file)
