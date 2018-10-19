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
from MaskedCCD import MaskedCCD
from BrightPixels import BrightPixels
from EOTestResults import EOTestResults
from generate_mask import generate_mask

from GlowingSources import GlowingSources

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
    nsig = pexConfig.Field("Number of sigma for bright source identification",
                           int, default=10)
    npix_min = pexConfig.Field("Minimum number of pixels for an extended bright source to be identified",
                               int, default=20)
    glowing_sources_file = pexConfig.Field("Glowing sources results filename",
                                           str, default='glowing_sources_params.fits')

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
        md = imutils.Metadata(dark_files[0], 1)
        exptime = ccd.md.get('EXPTIME')
        total_bright_pixels = 0
        total_bright_columns = 0
        if self.config.verbose:
            self.log.info("Amp         # bright pixels     # bright columns")
        #
        # Write bright pixel and column counts to results file.
        # Extended sources are found and saved to separate file.
        #
        results_file = self.config.eotest_results_file
        if results_file is None:
            results_file = os.path.join(self.config.output_dir,
                                        '%s_eotest_results.fits' % sensor_id)

        results = EOTestResults(results_file, namps=len(ccd))
        pixels = {}
        columns = {}
        nsig = self.config.nsig
        npix_min = self.config.npix_min
        outfile = self.config.glowing_sources_file
        glowing_sources = GlowingSources(nsig, npix_min)
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

            glowing_sources.find(ccd, amp)

        if self.config.verbose:
            self.log.info("Total bright pixels: %i" % total_bright_pixels)
            self.log.info("Total bright columns: %i" % total_bright_columns)
        glowing_sources.write_results(outfile)
        results.write(clobber=True)

        # Generate the mask file based on the pixel and columns.
        mask_file = os.path.join(self.config.output_dir,
                                 '%s_bright_pixel_mask.fits' % sensor_id)
        if os.path.isfile(mask_file):
            os.remove(mask_file)
        generate_mask(medfile, mask_file, self.config.mask_plane,
                      pixels=pixels, columns=columns)

if __name__ == '__main__':

    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument('sensor_id', type=str, 
                        help='Identifier for the sensor (e.g. S00, S01, etc)')
    parser.add_argument('dark_files', nargs='+',
                        help='Dark files to use for analysis.')
    parser.add_argument('--bias_frame', '-b', type=str, default=None)
    parser.add_argument('--eotest', '-e', type=str, default=None)
    args = parser.parse_args()

    sensor_id = args.sensor_id
    dark_files = args.dark_files
    mask_files = []
    if args.eotest is None:
        gains = dict((i, 1) for i in range(1, 17))
    bias_frame = args.bias_frame

    bright_pixels = BrightPixelsTask()
    bright_pixels.run(sensor_id, dark_files, tuple(mask_files), 
                      gains, bias_frame = bias_frame)
        
