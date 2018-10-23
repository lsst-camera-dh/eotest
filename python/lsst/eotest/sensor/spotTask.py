from future import print_function
from future import absolute_import
import os
import numpy as np
from astropy.io import fits

import lsst.eotest.image_utils as imutils
import lsst.afw.image as afwImage
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.meas.extensions.shapeHSM

from .MaskedCCD import MaskedCCD
from .AmplifierGeometry import parse_geom_kwd
from lsst.pipe.tasks.characterizeImage import CharacterizeImageTask, CharacterizeImageConfig


class SpotConfig(pexConfig.Config):
    """Configuration for Spot analysis task"""
    minpixels = pexConfig.Field("Minimum number of pixels above detection threshold", 
                                int, default=10)
    nsig = pexConfig.Field("Source footprint threshold in number of standard deviations of image section", 
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
    def run(self, sensor_id, infile, mask_files, bias_frame=None,
            spot_catalog=None, oscan_fit_order=1):

        imutils.check_temperatures(infiles, self.config.temp_set_point_tol,
                                   setpoint=self.config.temp_set_point,
                                   warn_only=True)

        if self.config.verbose and spot_catalog is None:
            self.log.info("Input files:")
            self.log.info("  {0}".format(infile))

        #
        # Detect and fit spots, accumulating the results by amplifier.
        #
        charConfig = CharacterizeImageConfig()
        charConfig.doMeasurePsf = False
        charConfig.doApCorr = False
        charConfig.repair.doCosmicRay = False
        charConfig.detection.minPixels = self.config.minpixels
        charConfig.detection.background.binSize = 10
        charConfig.detection.thresholdType = "stdev"
        charConfig.detection.thresholdValue = self.config.nsig
        charConfig.measurement.plugins.names |= ["ext_shapeHSM_HsmSourceMoments"]

        charTask = CharacterizeImageTask(config=charConfig)
        if self.config.verbose:
            self.log.info("processing {0}".format(infile))

        image = make_ccd_mosaic(infile)
        exposure = afwImage.ExposureF(image.getBBox())
        exposure.setImage(image)
                
        result = charTask.characterize(exposure)
        if self.config.verbose:
            self.log.info("Detected {0} objects".format(len(result.sourceCat)))

        if self.config.output_file is None:
            output_file = os.path.join(self.config.output_dir,
                                       '{0}_spot_results_nsig{1}.cat'.format(sensor_id, self.config.nsig))
        else:
            output_file = self.config.output_file
        if self.config.verbose:
            self.log.info("Writing spot results file to {0}".format(spot_results))
        src = result.sourceCat
        src.writeFits(output_file)
        
if __name__ == '__main__':

    import glob
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('sensor_id', help='Sensor name, e.g. S00')
    parser.add_argument('image_files', nargs='+', help='List of spot images')
    parser.add_argument('-b', '--bias_frame', default=None, 
                        help='Bias frame to use')
    parser.add_argument('-o', '--output_dir', 
                        default='/u/ec/elp25/private/spots_testing/',
                        help='Output directory')
    parser.add_argument('-n', '--nsig', type=float, default=10.0,
                        help='Number of standard deviations to use when setting threshold.')
    parser.add_argument('-m', '--mosaic', type=bool, default=False,
                        help='Boolean to mosaic CCD before spot finding.')
    args = parser.parse_args()

    sensor_id = args.sensor_id
    image_files = args.image_files
    bias_frame = args.bias_frame
    output_dir = args.output_dir
    nsig = args.nsig
    mosaic = args.mosaic
    mask_files = tuple()

    print image_files
    print bias_frame

    spottask = SpotTask()
    spottask.config.verbose = False
    spottask.config.output_dir = output_dir
    spottask.config.nsig = nsig
    spottask.run(sensor_id, image_files, mask_files, bias_frame=bias_frame,
                 mosaic=mosaic)
