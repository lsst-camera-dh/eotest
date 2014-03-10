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
    temp_tol = pexConfig.Field("Temperature tolerance in degrees C for CCDTEMP among dark files",
                               float, default=1.5)
    output_dir = pexConfig.Field("Output directory", str, default=".")
    eotest_results_file = pexConfig.Field("EO test results filename",
                                          str, default=None)
    verbose = pexConfig.Field("Turn verbosity on", bool, default=True)

class BrightPixelsTask(pipeBase.Task):
    """Task to find bright pixels and columns."""
    ConfigClass = BrightPixelsConfig
    _DefaultName = "BrightPixelsTask"

    @pipeBase.timeMethod
    def run(self, sensor_id, dark_files, mask_files, gains):
        imutils.check_temperatures(dark_files, self.config.temp_tol)
        median_images = {}
        for amp in imutils.allAmps:
            median_images[amp] = imutils.fits_median(dark_files,
                                                     imutils.dm_hdu(amp))
            medfile = os.path.join(self.config.output_dir,
                                   '%s_median_dark_bp.fits' % sensor_id)
        imutils.writeFits(median_images, medfile, dark_files[0])

        ccd = MaskedCCD(medfile, mask_files=mask_files)
        md = afwImage.readMetadata(dark_files[0], 1)
        exptime = ccd.md.get('EXPTIME')
        outfile = os.path.join(self.config.output_dir,
                               '%s_bright_pixel_map.fits' % sensor_id)
        total_bright_pixels = 0
        if self.config.verbose:
            self.log.info("Segment     # bright pixels")
        #
        # Write bright pixel counts to results file.
        #
        if self.config.eotest_results_file is None:
            results_file = '%s_eotest_results.fits' % sensor_id
        else:
            results_file = self.config.eotest_results_file
        
        results = EOTestResults(results_file)
        for amp in imutils.allAmps:
            bright_pixels = BrightPixels(ccd, amp, exptime, gains[amp])
            
            pixels, columns = bright_pixels.find()
            bright_pixels.generate_mask(outfile)
            count = len(pixels)
            total_bright_pixels += count
            results.add_seg_result(amp, 'NUM_BRIGHT_PIXELS', count)
            self.log.info("%s          %i" % (imutils.channelIds[amp], count))
        if self.config.verbose:
            self.log.info("Total bright pixels: %i" % total_bright_pixels)
        results.write(clobber=True)
