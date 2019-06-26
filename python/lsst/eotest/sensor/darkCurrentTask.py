"""
@brief Dark current task: compute 95th percentile dark current in
units of e-/sec/pixel.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
from __future__ import print_function
from __future__ import absolute_import
import os
import warnings
import psutil
import numpy as np
import astropy.io.fits as fits
from astropy.utils.exceptions import AstropyWarning, AstropyUserWarning
from lsst.eotest.fitsTools import fitsWriteto
import lsst.eotest.image_utils as imutils
from .MaskedCCD import MaskedCCD
from .EOTestResults import EOTestResults
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase


class DarkCurrentConfig(pexConfig.Config):
    """Configuration for DarkCurrentTask"""
    temp_set_point = pexConfig.Field("Required temperature (C) set point",
                                     float, default=-95.)
    temp_set_point_tol = pexConfig.Field("Required temperature set point tolerance (degrees C)",
                                         float, default=1.)
    output_dir = pexConfig.Field("Output directory", str, default=".")
    eotest_results_file = pexConfig.Field('EO test results filename',
                                          str, default=None)
    verbose = pexConfig.Field("Turn verbosity on", bool, default=True)


class DarkCurrentTask(pipeBase.Task):
    """Task to evaluate dark current quantiles."""
    ConfigClass = DarkCurrentConfig
    _DefaultName = "DarkCurrentTask"

    def log_mem_info(self, message, amp=None):
        if self.process is None:
            return
        rss_mem = self.process.memory_info().rss/1024.**3
        amp_info = ", amp %d" % amp if amp is not None else ''
        self.log.info('DarkCurrentTask, %s: rss_mem=%.2f GB; %s',
                      self.process.pid, rss_mem, message + amp_info)

    def _set_process(self):
        if os.environ.get('LCATR_LOG_MEM_INFO', False) == 'True':
            self.process = psutil.Process(os.getpid())
        else:
            self.process = None

    @pipeBase.timeMethod
    def run(self, sensor_id, dark_files, mask_files, gains, bias_frame=None):
        imutils.check_temperatures(dark_files, self.config.temp_set_point_tol,
                                   setpoint=self.config.temp_set_point,
                                   warn_only=True)
        median_images = {}
        self._set_process()
        self.log_mem_info("md = imutils.Metadata(...")
        md = imutils.Metadata(dark_files[0])
        for amp in imutils.allAmps(dark_files[0]):
            self.log_mem_info("imutils.fits_median", amp=amp)
            median_images[amp] = imutils.fits_median(dark_files,
                                                     imutils.dm_hdu(amp))
        medfile = os.path.join(self.config.output_dir,
                               '%s_median_dark_current.fits' % sensor_id)
        self.log_mem_info("imutils.writeFits(...")
        imutils.writeFits(median_images, medfile, dark_files[0])

        self.log_mem_info("ccd = MaskedCCD(...")
        ccd = MaskedCCD(medfile, mask_files=mask_files, bias_frame=bias_frame)

        dark95s = {}
        dark_medians = {}
        try:
            exptime = md.get('DARKTIME')
        except:
            exptime = md.get('EXPTIME')
        if self.config.verbose:
            self.log.info("Amp        95 percentile    median")
        dark_curr_pixels = []
        dark_curr_pixels_per_amp = {}
        for amp in ccd:
            imaging_region = ccd.amp_geom.imaging
            overscan = ccd.amp_geom.serial_overscan
            self.log_mem_info("image = imutils.unbias_and_trim", amp=amp)
            image = imutils.unbias_and_trim(im=ccd[amp].getImage(),
                                            overscan=overscan,
                                            imaging=imaging_region)
            self.log_mem_info("mask = imutils.trim", amp=amp)
            mask = imutils.trim(ccd[amp].getMask(), imaging_region)
            unmasked = []
            self.log_mem_info("for im_pixel, mask_pixel in...", amp=amp)
            for im_pixel, mask_pixel in zip(image.array.ravel(),
                                            mask.array.ravel()):
                if mask_pixel == 0:
                    unmasked.append(im_pixel)
            self.log_mem_info("unmasked.sort()", amp=amp)
            unmasked.sort()
            unmasked = np.array(unmasked)*gains[amp]/exptime
            dark_curr_pixels_per_amp[amp] = unmasked
            dark_curr_pixels.extend(unmasked)
            try:
                dark95s[amp] = unmasked[int(len(unmasked)*0.95)]
                dark_medians[amp] = unmasked[len(unmasked)//2]
            except IndexError as eobj:
                print(str(eobj))
                dark95s[amp] = -1.
                dark_medians[amp] = -1.
            if self.config.verbose:
                self.log.info("%2i         %.2e         %.2e"
                              % (amp, dark95s[amp], dark_medians[amp]))
        #
        # Compute 95th percentile dark current for CCD as a whole.
        #
        dark_curr_pixels = sorted(dark_curr_pixels)
        darkcurr95 = dark_curr_pixels[int(len(dark_curr_pixels)*0.95)]
        self.log_mem_info("dark95mean =")
        dark95mean = np.mean(list(dark95s.values()))
        if self.config.verbose:
            #self.log.info("CCD: mean 95 percentile value = %s" % dark95mean)
            self.log.info("CCD-wide 95 percentile value = %s" % darkcurr95)
        #
        # Update header of dark current median image file with dark
        # files used and dark95 values, and write dark95 values to the
        # eotest results file.
        #
        results_file = self.config.eotest_results_file
        if results_file is None:
            results_file = os.path.join(self.config.output_dir,
                                        '%s_eotest_results.fits' % sensor_id)
        results = EOTestResults(results_file, namps=len(ccd))
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=UserWarning, append=True)
            warnings.filterwarnings('ignore', category=AstropyWarning,
                                    append=True)
            warnings.filterwarnings('ignore', category=AstropyUserWarning,
                                    append=True)

            with fits.open(medfile) as output:
                for i, dark in enumerate(dark_files):
                    output[0].header['DARK%02i' % i] = os.path.basename(dark)
                # Write overall dark current 95th percentile
                results.output['AMPLIFIER_RESULTS'].header['DARK95'] \
                    = darkcurr95
                for amp in ccd:
                    output[0].header['DARK95%s'%imutils.channelIds[amp]] \
                        = dark95s[amp]
                    results.add_seg_result(amp, 'DARK_CURRENT_95', dark95s[amp])
                    results.add_seg_result(amp, 'DARK_CURRENT_MEDIAN',
                                           dark_medians[amp])
                fitsWriteto(output, medfile, overwrite=True, checksum=True)
                results.write(clobber=True)

        return dark_curr_pixels_per_amp, dark95s
