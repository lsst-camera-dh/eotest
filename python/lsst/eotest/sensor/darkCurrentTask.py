"""
@brief Dark current task: compute 95th percentile dark current in
units of e-/sec/pixel.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
from collections import defaultdict
import warnings
import numpy as np
import pandas as pd
import astropy.io.fits as fits
from astropy.utils.exceptions import AstropyWarning, AstropyUserWarning
import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.utils.timer import timeMethod
from lsst.eotest.fitsTools import fitsWriteto
import lsst.eotest.image_utils as imutils
from .MaskedCCD import MaskedCCD
from .EOTestResults import EOTestResults


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

    @timeMethod
    def run(self, sensor_id, dark_files, mask_files, gains, bias_frame=None,
            linearity_correction=None, dark_files_linear_fit=None):
        if dark_files_linear_fit is not None:
            self.linear_fit(sensor_id, dark_files_linear_fit, mask_files,
                            gains, bias_frame=bias_frame,
                            linearity_correction=linearity_correction)
        return self.compute_percentiles(sensor_id, dark_files, mask_files,
                                        gains, bias_frame=bias_frame,
                                        linearity_correction=linearity_correction)

    def linear_fit(self, sensor_id, dark_files, mask_files, gains,
                   bias_frame=None, linearity_correction=None):
        """Fit a slope and intercept for a dark frame data set that has
        a range of different exposure times.
        """
        data = defaultdict(list)
        for item in dark_files:
            ccd = MaskedCCD(item, mask_files=mask_files, bias_frame=bias_frame,
                            linearity_correction=linearity_correction)
            for amp in ccd:
                data['amp'].append(amp)
                data['darktime'].append(ccd.md.get('DARKTIME'))
                image = ccd.unbiased_and_trimmed_image(amp)
                stats = afwMath.makeStatistics(image, afwMath.MEDIAN)
                data['median'].append(stats.getValue(afwMath.MEDIAN)*gains[amp])
        df0 = pd.DataFrame(data=data)
        slopes, intercepts = dict(), dict()
        for amp in ccd:
            df = df0.query(f'amp == {amp}')
            slopes[amp], intercepts[amp] \
                = np.polyfit(df['darktime'], df['median'], 1)

        results_file = self.config.eotest_results_file
        if results_file is None:
            results_file = os.path.join(self.config.output_dir,
                                        '%s_eotest_results.fits' % sensor_id)
        results = EOTestResults(results_file, namps=len(ccd))
        # Write slopes and intercepts for each amp
        for amp in ccd:
            results.add_seg_result(amp, 'DARK_CURRENT_SLOPE', slopes[amp])
            results.add_seg_result(amp, 'DARK_CURRENT_INTERCEPT',
                                   intercepts[amp])
        results.write(clobber=True)

    def compute_percentiles(self, sensor_id, dark_files, mask_files, gains,
                            bias_frame=None, linearity_correction=None):
        """Compute median and 95th percentiles of pixel values for dark frame
        data sets that have the same integration times."""
        imutils.check_temperatures(dark_files, self.config.temp_set_point_tol,
                                   setpoint=self.config.temp_set_point,
                                   warn_only=True)
        median_images = {}
        md = imutils.Metadata(dark_files[0])
        for amp in imutils.allAmps(dark_files[0]):
            median_images[amp] = imutils.fits_median(dark_files,
                                                     imutils.dm_hdu(amp))
        medfile = os.path.join(self.config.output_dir,
                               '%s_median_dark_current.fits' % sensor_id)
        imutils.writeFits(median_images, medfile, dark_files[0])

        ccd = MaskedCCD(medfile, mask_files=mask_files, bias_frame=bias_frame,
                        linearity_correction=linearity_correction)

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
            mi = ccd.unbiased_and_trimmed_image(amp)
            imarr = mi.getImage().array
            mskarr = mi.getMask().array
            pixels = imarr.reshape(1, imarr.shape[0]*imarr.shape[1])[0]
            masked = mskarr.reshape(1, mskarr.shape[0]*mskarr.shape[1])[0]
            unmasked = [pixels[i] for i in range(len(pixels)) if masked[i] == 0]
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
