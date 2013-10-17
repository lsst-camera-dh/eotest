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
from fe55_gain_fitter import fe55_gain_fitter

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

class Fe55Config(pexConfig.Config):
    """Configuration for Fe55 analysis task"""
    chiprob_min = pexConfig.Field("Minimum chi-square probability for cluster fit",
                                  float, default=0.1)
    output_dir = pexConfig.Field("Output directory", str, default='.')
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
        fitter = PsfGaussFit()
        for infile in infiles:
            if self.config.verbose:
                self.log.info("processing %s" % infile)
            ccd = MaskedCCD(infile, mask_files=mask_files)
            for amp in ccd:
                if self.config.verbose:
                    self.log.info("  amp %i" % amp)
                fitter.process_image(ccd, amp)
        psf_results = os.path.join(self.config.output_dir,
                                   '%s_psf_results.fits' % sensor_id)
        fitter.write_results(outfile=psf_results)
        #
        # Fit the DN distributions to obtain the system gain per amp.
        #
        # @todo Find a better way of saving the gain estimates.
        #
        gains = {}
        output = pyfits.HDUList()
        output.append(pyfits.PrimaryHDU())
        for amp in imutils.allAmps:
            data = fitter.results(min_prob=self.config.chiprob_min, amp=amp)
            dn = data[1]
            gains[amp] = fe55_gain_fitter(dn, make_plot=False)
            output[0].header.update('GAIN%s' % imutils.channelIds[amp], gains[amp])

        gain_file = os.path.join(self.config.output_dir, '%s_gains.fits' % sensor_id)
        output.writeto(gain_file, clobber=True, checksum=True)

        return gains
