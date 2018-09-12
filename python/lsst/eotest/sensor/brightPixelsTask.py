"""
@brief Bright pixels task: Find pixels and columns in a median image
constructed from an ensemble of darks.  The bright pixel threshold is
specified via the --ethresh option and is in units of -e per pixel per
second.  The threshold for the number of bright pixels that define a
bright column is specified via the --colthresh option.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os

import lsst.eotest.image_utils as imutils
from .MaskedCCD import MaskedCCD
from .BrightPixels import BrightPixels
from .EOTestResults import EOTestResults
from .generate_mask import generate_mask

import lsst.afw.image as afwImage
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase


class BrightPixelsConfig(pexConfig.Config):
    """Configuration for bright pixels task"""
    ethresh = pexConfig.Field("Bright pixel threshold in e- per pixel per second",
                              int, default=5)
    colthresh = pexConfig.Field("Bright column threshold in # bright pixels",
                                int, default=20)
    mask_plane = pexConfig.Field("Mask plane to be used for output mask file",
                                 str, default='BAD')
    temp_set_point = pexConfig.Field("Required temperature (C) set point",
                                     float, default=-95.)
    temp_set_point_tol = pexConfig.Field("Required temperature set point tolerance (degrees C)",
                                         float, default=1.)
    output_dir = pexConfig.Field("Output directory", str, default=".")
    eotest_results_file = pexConfig.Field("EO test results filename",
                                          str, default=None)
    verbose = pexConfig.Field("Turn verbosity on", bool, default=True)


class BrightPixelsTask(pipeBase.Task):
    """Task to find bright pixels and columns."""
    ConfigClass = BrightPixelsConfig
    _DefaultName = "BrightPixelsTask"

    @pipeBase.timeMethod
    def run(self, sensor_id, dark_files, mask_files, gains, bias_frame=None):
        imutils.check_temperatures(dark_files, self.config.temp_set_point_tol,
                                   setpoint=self.config.temp_set_point,
                                   warn_only=True)
        median_images = {}
        for amp in imutils.allAmps(dark_files[0]):
            median_images[amp] = imutils.fits_median(dark_files,
                                                     imutils.dm_hdu(amp))
        medfile = os.path.join(self.config.output_dir,
                               '%s_median_dark_bp.fits' % sensor_id)
        imutils.writeFits(median_images, medfile, dark_files[0])

        ccd = MaskedCCD(medfile, mask_files=mask_files, bias_frame=bias_frame)
        md = imutils.Metadata(dark_files[0])
        exptime = ccd.md.get('EXPTIME')
        total_bright_pixels = 0
        total_bright_columns = 0
        if self.config.verbose:
            self.log.info("Amp         # bright pixels     # bright columns")
        #
        # Write bright pixel and column counts to results file.
        #
        results_file = self.config.eotest_results_file
        if results_file is None:
            results_file = os.path.join(self.config.output_dir,
                                        '%s_eotest_results.fits' % sensor_id)

        results = EOTestResults(results_file, namps=len(ccd))
        pixels = {}
        columns = {}
        for amp in ccd:
            bright_pixels = BrightPixels(ccd, amp, exptime, gains[amp])
            pixels[amp], columns[amp] = bright_pixels.find()
            pix_count = len(pixels[amp])
            col_count = len(columns[amp])
            total_bright_pixels += pix_count
            total_bright_columns += col_count
            results.add_seg_result(amp, 'NUM_BRIGHT_PIXELS', pix_count)
            results.add_seg_result(amp, 'NUM_BRIGHT_COLUMNS', col_count)
            self.log.info("%2i          %i          %i" %
                          (amp, pix_count, col_count))
        if self.config.verbose:
            self.log.info("Total bright pixels: %i" % total_bright_pixels)
            self.log.info("Total bright columns: %i" % total_bright_columns)
        results.write(overwrite=True)

        # Generate the mask file based on the pixel and columns.
        mask_file = os.path.join(self.config.output_dir,
                                 '%s_bright_pixel_mask.fits' % sensor_id)
        if os.path.isfile(mask_file):
            os.remove(mask_file)
        generate_mask(medfile, mask_file, self.config.mask_plane,
                      pixels=pixels, columns=columns)
