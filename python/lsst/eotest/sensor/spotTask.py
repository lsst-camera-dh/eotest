import os
import numpy as np
from astropy.io import fits
import lsst.eotest.image_utils as imutils
from spot_psf import SpotMomentFit
from MaskedCCD import MaskedCCD
from EOTestResults import EOTestResults

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

class SpotConfig(pexConfig.Config):
    """Configuration for Spot analysis task"""
    chiprob_min = pexConfig.Field("Minimum chi-square probability for cluster fit",
                                  float, default=0.1)
    nsig = pexConfig.Field("Charge cluster footprint threshold in number of standard deviations of noise in bias section", float, default=2)
    temp_set_point = pexConfig.Field("Required temperature (C) set point",
                                     float, default=-95.)
    temp_set_point_tol = pexConfig.Field("Required temperature set point tolerance (degrees C)",
                                         float, default=1.)

    output_dir = pexConfig.Field("Output directory", str, default='.')
    output_file = pexConfig.Field("Output filename", str, default=None)
    eotest_results_file = pexConfig.Field("EO test results filename",
                                          str, default=None)
    verbose = pexConfig.Field("Turn verbosity on", bool, default=True)

class SpotTask(pipeBase.Task):
    """Task to estimate spot moments from spot projector data."""

    ConfigClass = SpotConfig
    _DefaultName = "SpotTask"

    @pipeBase.timeMethod
    def run(self, sensor_id, infiles, mask_files, bias_frame=None,
            spot_catalog=None, chiprob_min=0.1, accuracy_req=0,
            oscan_fit_order=1):

        imutils.check_temperatures(infiles, self.config.temp_set_point_tol,
                                   setpoint=self.config.temp_set_point,
                                   warn_only=True)

        if self.config.verbose and spot_catalog is None:
            self.log.info("Input files:")
            for item in infiles:
                self.log.info("  %s" % item)

        #
        # Detect and fit spots, accumulating the results by amplifier.
        #
        fitter = SpotMomentFit(nsig=self.config.nsig)
        if spot_catalog is None:
            for infile in infiles:
                if self.config.verbose:
                    self.log.info("processing %s" % infile)

                ccd = MaskedCCD(infile, mask_files=mask_files,
                                bias_frame=bias_frame)
                for amp in ccd:
                    if self.config.verbose:
                        self.log.info("  amp %i" % amp)
                    fitter.process_image(ccd, amp, oscan_fit_order=oscan_fit_order)

            if self.config.output_file is None:
                spot_results = os.path.join(self.config.output_dir,
                                           '%s_spot_results_nsig%i.fits'
                                           % (sensor_id, self.config.nsig))
            else:
                spot_results = self.config.output_file
            if self.config.verbose:
                self.log.info("Writing spot results file to %s" % spot_results)
            fitter.write_results(outfile=spot_results)
            namps = len(ccd)
        else:
            fitter.read_spot_catalog(spot_catalog)
            namps = fits.open(spot_catalog)[0].header['NAMPS']

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
    args = parser.parse_args()

    sensor_id = args.sensor_id
    image_files = args.image_files
    bias_frame = args.bias_frame
    output_dir = args.output_dir
    nsig = args.nsig
    mask_files = tuple()

    print image_files
    print bias_frame

    spottask = SpotTask()
    spottask.config.verbose = False
    spottask.config.output_dir = output_dir
    spottask.config.nsig = nsig
    spottask.run(sensor_id, image_files, mask_files, bias_frame)
