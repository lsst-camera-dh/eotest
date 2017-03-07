"""
@brief PSF characterization and system gain from distribution of
Gaussian fit parameters to Fe55 data.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import numpy as np
import astropy.io.fits as fits
import lsst.eotest.image_utils as imutils
from fe55_psf import PsfGaussFit, psf_sigma_statistics
from MaskedCCD import MaskedCCD
from EOTestResults import EOTestResults
from Fe55GainFitter import Fe55GainFitter

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

class Fe55Config(pexConfig.Config):
    """Configuration for Fe55 analysis task"""
    chiprob_min = pexConfig.Field("Minimum chi-square probability for cluster fit",
                                  float, default=0.1)
    nsig = pexConfig.Field("Charge cluster footprint threshold in number of standard deviations of noise in bias section", float, default=4)
    temp_set_point = pexConfig.Field("Required temperature (C) set point",
                                     float, default=-95.)
    temp_set_point_tol = pexConfig.Field("Required temperature set point tolerance (degrees C)",
                                         float, default=1.)

    output_dir = pexConfig.Field("Output directory", str, default='.')
    output_file = pexConfig.Field("Output filename", str, default=None)
    eotest_results_file = pexConfig.Field("EO test results filename",
                                          str, default=None)
    fit_xy = pexConfig.Field("Fit sigmas in x- and y-directions separately",
                             bool, default=False)
    verbose = pexConfig.Field("Turn verbosity on", bool, default=True)

class Fe55Task(pipeBase.Task):
    """Task to estimate PSF size and system gain from the distribution of
    Gaussian fit parameters to Fe55 data."""
    ConfigClass = Fe55Config
    _DefaultName = "Fe55Task"

    def fit_gains(self, fitter, gains, gain_errors, sigma_modes, amps=None):
        "Fit the DN distributions to obtain the system gain per amp."
        my_gains, my_gain_errors, my_sigma_modes = \
            gains, gain_errors, sigma_modes
        if amps is None:
            amps = imutils.allAmps()
        for amp in amps:
            data = fitter.results(min_prob=self.config.chiprob_min, amp=amp)
            dn = data['dn']
            if len(dn) > 2:
                try:
                    foo = Fe55GainFitter(dn)
                    kalpha_peak, kalpha_sigma = foo.fit()
                    my_gains[amp] = foo.gain
                    my_gain_errors[amp] = foo.gain_error
                except RuntimeError as eobj:
                    print(eobj)
                    continue
                try:
                    sigma = sorted(np.concatenate((data['sigmax'],
                                                   data['sigmay']))*10)
                    mode, median, mean = psf_sigma_statistics(sigma, bins=50,
                                                              range=(2, 6))
                    my_sigma_modes[amp] = float(mode)
                except RuntimeError as eobj:
                    print(eobj)
                    continue
        return my_gains, my_gain_errors, my_sigma_modes

    @pipeBase.timeMethod
    def run(self, sensor_id, infiles, mask_files, bias_frame=None,
            fe55_catalog=None, minClustersPerAmp=None, chiprob_min=0.1,
            accuracy_req=0, oscan_fit_order=1):
        imutils.check_temperatures(infiles, self.config.temp_set_point_tol,
                                   setpoint=self.config.temp_set_point,
                                   warn_only=True)
        if self.config.verbose and fe55_catalog is None:
            self.log.info("Input files:")
            for item in infiles:
                self.log.info("  %s" % item)
        #
        # Detect and fit 2D Gaussian to Fe55 charge clusters,
        # accumulating the results by amplifier.
        #
        fitter = PsfGaussFit(nsig=self.config.nsig, fit_xy=self.config.fit_xy)
        gains, gain_errors, sigma_modes = {}, {}, {}
        if fe55_catalog is None:
            for infile in infiles:
                if self.config.verbose:
                    self.log.info("processing %s" % infile)
                ccd = MaskedCCD(infile, mask_files=mask_files,
                                bias_frame=bias_frame)
                for amp in ccd:
                    if self.config.verbose:
                        self.log.info("  amp %i" % amp)
                    if gain_errors.has_key(amp):
                        gain_accuracy = np.abs(gain_errors[amp]/gains[amp])
                        if self.config.verbose:
                            message = "  Relative gain accuracy, dgain/gain " \
                                 + "= %.2e" % gain_accuracy
                            self.log.info(message)
                        if gain_accuracy < accuracy_req:
                            # Requested accuracy already obtained, so
                            # skip cluster fitting.
                            continue
                    fitter.process_image(ccd, amp, logger=self.log,
                                         oscan_fit_order=oscan_fit_order)
                    gains, gain_errors, sigma_modes = \
                        self.fit_gains(fitter, gains, gain_errors, sigma_modes,
                                       amps=ccd.keys())
            if self.config.output_file is None:
                psf_results = os.path.join(self.config.output_dir,
                                           '%s_psf_results_nsig%i.fits'
                                           % (sensor_id, self.config.nsig))
            else:
                psf_results = self.config.output_file
            if self.config.verbose:
                self.log.info("Writing psf results file to %s" % psf_results)
            fitter.write_results(outfile=psf_results)
            namps = len(ccd)
        else:
            fitter.read_fe55_catalog(fe55_catalog)
            namps = fits.open(fe55_catalog)[0].header['NAMPS']

        results_file = self.config.eotest_results_file
        if results_file is None:
            results_file = os.path.join(self.config.output_dir,
                                        '%s_eotest_results.fits' % sensor_id)

        if self.config.verbose:
            self.log.info("Writing gain and psf sigma results to %s"
                          % results_file)

        results = EOTestResults(results_file, namps=namps)
        for amp in gains:
            results.add_seg_result(amp, 'GAIN', gains[amp])
            results.add_seg_result(amp, 'GAIN_ERROR', gain_errors[amp])
            results.add_seg_result(amp, 'PSF_SIGMA', sigma_modes[amp])
        results.write(clobber=True)

        return gains
