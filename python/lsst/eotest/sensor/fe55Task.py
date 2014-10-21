"""
@brief PSF characterization and system gain from distribution of
Gaussian fit parameters to Fe55 data.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import numpy as np
import pyfits
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

    @pipeBase.timeMethod
    def run(self, sensor_id, infiles, mask_files, bias_frame=None,
            fe55_catalog=None):
        if self.config.verbose and fe55_catalog is None:
            self.log.info("Input files:")
            for item in infiles:
                self.log.info("  %s" % item) 
        #
        # Detect and fit 2D Gaussian to Fe55 charge clusters,
        # accumulating the results by amplifier.
        #
        fitter = PsfGaussFit(nsig=self.config.nsig, fit_xy=self.config.fit_xy)
        if fe55_catalog is None:
            for infile in infiles:
                if self.config.verbose:
                    self.log.info("processing %s" % infile)
                ccd = MaskedCCD(infile, mask_files=mask_files,
                                bias_frame=bias_frame)
                for amp in ccd:
                    if self.config.verbose:
                        self.log.info("  amp %i" % amp)
                    fitter.process_image(ccd, amp, logger=self.log)
            if self.config.output_file is None:
                psf_results = os.path.join(self.config.output_dir,
                                           '%s_psf_results_nsig%i.fits' 
                                           % (sensor_id, self.config.nsig))
            else:
                psf_results = self.config.output_file
            if self.config.verbose:
                self.log.info("Writing psf results file to %s" % psf_results)
            fitter.write_results(outfile=psf_results)
        else:
            fitter.read_fe55_catalog(fe55_catalog)
        #
        # Fit the DN distributions to obtain the system gain per amp.
        #
        gains = {}
        sigma_modes = {}
        for amp in imutils.allAmps:
            data = fitter.results(min_prob=self.config.chiprob_min, amp=amp)
            dn = data['dn']
            if len(dn) > 2:
                try:
                    foo = Fe55GainFitter(dn)
                    kalpha_peak, kalpha_sigma = foo.fit()
                    gains[amp] = foo.gain
                except RuntimeError, e:
                    print e
                    continue
                try:
                    sigma = np.concatenate((data['sigmax'], data['sigmay']))*10
                    mode, median, mean = psf_sigma_statistics(sigma, bins=50,
                                                              range=(2,6))
                    sigma_modes[amp] = float(mode)
                except RuntimeError, e:
                    print e
                    continue
            else:
                if self.config.verbose:
                    self.log.info("Too few charge clusters (%i) found for amp %s" % (len(dn), amp))

        results_file = self.config.eotest_results_file
        if results_file is None:
            results_file = os.path.join(self.config.output_dir, 
                                        '%s_eotest_results.fits' % sensor_id)
            
        if self.config.verbose:
            self.log.info("Writing gain and psf sigma results to %s" 
                          % results_file)

        results = EOTestResults(results_file)
        for amp in gains:
            results.add_seg_result(amp, 'GAIN', gains[amp])
            results.add_seg_result(amp, 'PSF_SIGMA', sigma_modes[amp])
        results.write(clobber=True)

        return gains
