"""
@brief PSF characterization and system gain from distribution of
Gaussian fit parameters to Fe55 data.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import numpy as np
import pyfits
import lsst.eotest.image_utils as imutils
from fe55_psf import PsfGaussFit
from MaskedCCD import MaskedCCD
from EOTestResults import EOTestResults
from fe55_gain_fitter import fe55_gain_fitter

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

class Fe55Config(pexConfig.Config):
    """Configuration for Fe55 analysis task"""
    chiprob_min = pexConfig.Field("Minimum chi-square probability for cluster fit",
                                  float, default=0.1)
    nsig = pexConfig.Field("Charge cluster footprint threshold in number of standard deviations of noise in bias section", float, default=4)
    output_dir = pexConfig.Field("Output directory", str, default='.')
    output_file = pexConfig.Field("Output filename", str, default=None)
    verbose = pexConfig.Field("Turn verbosity on", bool, default=True)
    
class Fe55Task(pipeBase.Task):
    """Task to estimate PSF size and system gain from the distribution of
    Gaussian fit parameters to Fe55 data."""
    ConfigClass = Fe55Config
    _DefaultName = "Fe55Task"

    @pipeBase.timeMethod
    def run(self, sensor_id, infiles, mask_files):
        #
        # Detect and fit 2D Gaussian to Fe55 charge clusters,
        # accumulating the results by amplifier.
        #
        fitter = PsfGaussFit(nsig=self.config.nsig)
        for infile in infiles:
            if self.config.verbose:
                self.log.info("processing %s" % infile)
            ccd = MaskedCCD(infile, mask_files=mask_files)
            for amp in ccd:
                if self.config.verbose:
                    self.log.info("  amp %i" % amp)
                fitter.process_image(ccd, amp)
        if self.config.output_file is None:
            psf_results = os.path.join(self.config.output_dir,
                                       '%s_psf_results.fits' % sensor_id)
        else:
            psf_results = self.config.output_file
        if self.config.verbose:
            self.log.info("Writing psf results file to %s" % psf_results)
        fitter.write_results(outfile=psf_results)
        #
        # Fit the DN distributions to obtain the system gain per amp.
        #
        gains = {}
        for amp in imutils.allAmps:
            data = fitter.results(min_prob=self.config.chiprob_min, amp=amp)
            dn = data[1]
            gains[amp], kalpha_peak, kalpha_sigma = fe55_gain_fitter(dn)
        results_file = os.path.join(self.config.output_dir, 
                                    '%s_results.fits' % sensor_id)
        if self.config.verbose:
            self.log.info("Writing gain file to %s" % results_file)

        results = EOTestResults(results_file)
        for amp in gains:
            results.add_seg_result(amp, 'GAIN', gains[amp])
        results.write(clobber=True)

        return gains
