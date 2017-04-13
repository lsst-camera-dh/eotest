"""
@brief Dark current task: compute 95th percentile dark current in
units of e-/sec/pixel.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import numpy as np
import astropy.io.fits as fits
from lsst.eotest.fitsTools import fitsWriteto
import lsst.eotest.image_utils as imutils
from MaskedCCD import MaskedCCD
from EOTestResults import EOTestResults
import lsst.afw.image as afwImage
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

    @pipeBase.timeMethod
    def run(self, sensor_id, dark_files, mask_files, gains, bias_frame=None):
        imutils.check_temperatures(dark_files, self.config.temp_set_point_tol,
                                   setpoint=self.config.temp_set_point,
                                   warn_only=True)
        median_images = {}
        md = imutils.Metadata(dark_files[0], 1)
        for amp in imutils.allAmps(dark_files[0]):
            median_images[amp] = imutils.fits_median(dark_files,
                                                     imutils.dm_hdu(amp))
        medfile = os.path.join(self.config.output_dir,
                               '%s_median_dark_current.fits' % sensor_id)
        imutils.writeFits(median_images, medfile, dark_files[0])

        ccd = MaskedCCD(medfile, mask_files=mask_files, bias_frame=bias_frame)

        dark95s = {}
        exptime = md.get('EXPTIME')
        if self.config.verbose:
            self.log.info("Amp        95 percentile    median")
        dark_curr_pixels = []
        dark_curr_pixels_per_amp = {}
        for amp in ccd:
            imaging_region = ccd.amp_geom.imaging
            overscan = ccd.amp_geom.serial_overscan
            image = imutils.unbias_and_trim(ccd[amp].getImage(),
                                            overscan, imaging_region)
            mask = imutils.trim(ccd[amp].getMask(), imaging_region)
            imarr = image.getArray()
            mskarr = mask.getArray()
            pixels = imarr.reshape(1, imarr.shape[0]*imarr.shape[1])[0]
            masked = mskarr.reshape(1, mskarr.shape[0]*mskarr.shape[1])[0]
            unmasked = [pixels[i] for i in range(len(pixels)) if masked[i] == 0]
            unmasked.sort()
            unmasked = np.array(unmasked)*gains[amp]/exptime
            dark_curr_pixels_per_amp[amp] = unmasked
            dark_curr_pixels.extend(unmasked)
            try:
                dark95s[amp] = unmasked[int(len(unmasked)*0.95)]
                median = unmasked[len(unmasked)/2]
            except IndexError as eobj:
                print str(eobj)
                dark95s[amp] = -1.
                median = -1.
            if self.config.verbose:
                self.log.info("%2i         %.2e         %.2e"
                              % (amp, dark95s[amp], median))
        #
        # Compute 95th percentile dark current for CCD as a whole.
        #
        dark_curr_pixels = sorted(dark_curr_pixels)
        darkcurr95 = dark_curr_pixels[int(len(dark_curr_pixels)*0.95)]
        dark95mean = np.mean(dark95s.values())
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
        output = fits.open(medfile)
        for i, dark in enumerate(dark_files):
            output[0].header['DARK%02i' % i] = os.path.basename(dark)
        # Write overall dark current 95th percentile
        results.output['AMPLIFIER_RESULTS'].header['DARK95'] = darkcurr95
        for amp in ccd:
            output[0].header['DARK95%s'%imutils.channelIds[amp]] = dark95s[amp]
            results.add_seg_result(amp, 'DARK_CURRENT_95', dark95s[amp])
        fitsWriteto(output, medfile, clobber=True, checksum=True)
        results.write(clobber=True)
        return dark_curr_pixels_per_amp, dark95s
