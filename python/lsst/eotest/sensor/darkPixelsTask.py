"""
@brief DarkPixelsTask: Using a medianed superflat dataset (at 500nm),
find the pixels that have less than 80% (or some other specified
fraction) of the mean (*not median*) signal among all the pxiels in a
given segment, excluding pixels masked as bright pixels or in the edge
rolloff/blooming stop areas.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os

import lsst.eotest.image_utils as imutils
from MaskedCCD import MaskedCCD
from DarkPixels import DarkPixels
from EOTestResults import EOTestResults
from cteTask import superflat

import lsst.afw.image as afwImage
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

class DarkPixelsConfig(pexConfig.Config):
    """Configuration for bright pixels task"""
    thresh = pexConfig.Field("Fractional threshold of segment mean for identifying dark pixels",
                             float, default=0.8)
    colthresh = pexConfig.Field("Threshold number of continguous dark pixels that a dark column",
                                int, default=100)
    mask_plane = pexConfig.Field("Mask plane to be used for output mask file",
                                 str, default='BAD')
    output_dir = pexConfig.Field("Output directory", str, default=".")
    eotest_results_file = pexConfig.Field("EO test results filename",
                                          str, default=None)
    verbose = pexConfig.Field("Turn verbosity on", bool, default=True)

class DarkPixelsTask(pipeBase.Task):
    """Task to find dark pixels and columns."""
    ConfigClass = DarkPixelsConfig
    _DefaultName = "DarkPixelsTask"

    @pipeBase.timeMethod
    def run(self, sensor_id, sflat_files, mask_files, bias_frame=None):
        medfile = os.path.join(self.config.output_dir,
                               '%s_median_sflat.fits' % sensor_id)
        superflat(sflat_files, outfile=medfile)

        ccd = MaskedCCD(medfile, mask_files=mask_files, bias_frame=bias_frame)
        md = imutils.Metadata(sflat_files[0], 1)
        outfile = os.path.join(self.config.output_dir,
                               '%s_dark_pixel_mask.fits' % sensor_id)
        if os.path.isfile(outfile):
            os.remove(outfile)
        total_dark_pixels = 0
        total_dark_columns = 0
        if self.config.verbose:
            self.log.info("Amp         # dark pixels     # dark columns")
        #
        # Write dark pixel and column counts to results file.
        #
        results_file = self.config.eotest_results_file
        if results_file is None:
            results_file = os.path.join(self.config.output_dir,
                                        '%s_eotest_results.fits' % sensor_id)

        results = EOTestResults(results_file)
        for amp in ccd:
            dark_pixels = DarkPixels(ccd, amp,
                                     frac_thresh=self.config.thresh,
                                     colthresh=self.config.colthresh,
                                     mask_plane=self.config.mask_plane)
            pixels, columns = dark_pixels.find()
            dark_pixels.generate_mask(outfile)
            pix_count = len(pixels)
            col_count = len(columns)
            total_dark_pixels += pix_count
            total_dark_columns += col_count
            results.add_seg_result(amp, 'NUM_DARK_PIXELS', pix_count)
            results.add_seg_result(amp, 'NUM_DARK_COLUMNS', col_count)
            self.log.info("%2i          %i          %i" 
                          % (amp, pix_count, col_count))
        if self.config.verbose:
            self.log.info("Total dark pixels: %i" % total_dark_pixels)
            self.log.info("Total dark columns: %i" % total_dark_columns)
        results.write(clobber=True)
