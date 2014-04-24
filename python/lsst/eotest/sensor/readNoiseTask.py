"""
@brief Compute read noise distributions for a sample of images.  Bias
and system readout noise exposures, the latter for determining the
noise contribution from the electronics, must be provided.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import numpy as np
import pyfits
import lsst.eotest.image_utils as imutils
from MaskedCCD import MaskedCCD
from EOTestResults import EOTestResults
from read_noise import noise_dists

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

def _write_read_noise_dists(outfile, Ntot, Nsys, gains, bias, sysnoise):
    output = pyfits.HDUList()
    output.append(pyfits.PrimaryHDU())
    output[0].header.update('BIASFILE', bias)
    output[0].header.update('SYSNFILE', sysnoise)
    for amp in imutils.allAmps:
        sigtot, sigsys = Ntot[amp], Nsys[amp]
        nread_col = pyfits.Column(name="TOTAL_NOISE", format="E",
                                  unit="e- rms", array=sigtot)
        nsys_col = pyfits.Column(name="SYSTEM_NOISE", format="E",
                                 unit="e- rms", array=sigsys)
        output.append(pyfits.new_table((nread_col, nsys_col)))
        output[amp].name = "AMP%s" % imutils.channelIds[amp]
        output[0].header.update("GAIN%s" % imutils.channelIds[amp], gains[amp])
        output[0].header.update("SIGTOT%s" % imutils.channelIds[amp],
                                imutils.median(sigtot))
        output[0].header.update("SIGSYS%s" % imutils.channelIds[amp],
                                imutils.median(sigsys))
    output.writeto(outfile, clobber=True)

class ReadNoiseConfig(pexConfig.Config):
    """Configuration for read noise task"""
    dx = pexConfig.Field("Subregion size in pixels along x-direction",
                         int, default=100)
    dy = pexConfig.Field("Subregion size in pixels along y-direction",
                         int, default=100)
    nsamp = pexConfig.Field("Number of subregions to sample",
                            int, default=1000)
    output_dir = pexConfig.Field("Output directory", str, default=".")
    eotest_results_file = pexConfig.Field("EO test results filename",
                                          str, default=None)
    verbose = pexConfig.Field("Turn verbosity on", bool, default=True)

class ReadNoiseTask(pipeBase.Task):
    """Task to estimate sensor read noise."""
    ConfigClass = ReadNoiseConfig
    _DefaultName = "ReadNoiseTask"

    @pipeBase.timeMethod
    def run(self, sensor_id, bias_files, system_noise_files, mask_files, gains):
        outfiles = []
        Nread_dists = dict([(amp, []) for amp in imutils.allAmps])
        for i, bias, sysnoise in zip(range(len(bias_files)), bias_files,
                                     system_noise_files):
            outfile = "%s_read_noise_%03i.fits" % (sensor_id, i)
            outfile = os.path.join(self.config.output_dir, outfile)
            outfiles.append(outfile)
            
            if self.config.verbose:
                self.log.info("Processing %s %s -> %s" % (bias, sysnoise,
                                                          outfile))

            # Determine the nominal imaging region from the bias file.
            ccd = MaskedCCD(bias)
            imaging = ccd.amp_geom.imaging
            #
            # Create a single sub-region sampler so that the same
            # sub-regions will be used for both the bias and system
            # noise frames.
            #
            sampler = imutils.SubRegionSampler(self.config.dx, self.config.dy,
                                               self.config.nsamp,
                                               imaging=imaging)

            Ntot = noise_dists(bias, gains, sampler, mask_files=mask_files)
            Nsys = noise_dists(sysnoise, gains, sampler, mask_files=mask_files)

            _write_read_noise_dists(outfile, Ntot, Nsys, gains, bias, sysnoise)
            
        results_file = self.config.eotest_results_file
        if results_file is None:
            results_file = os.path.join(self.config.output_dir,
                                        '%s_eotest_results.fits' % sensor_id)

        results = EOTestResults(results_file)
        if self.config.verbose:
            self.log.info("Amp    read noise")
        for amp in imutils.allAmps:
            var = imutils.median(Ntot[amp])**2 - imutils.median(Nsys[amp])**2
            if var >= 0:
                Nread = np.sqrt(var)
            else:
                Nread = -1
            line = "%s        %.4f" % (amp, Nread)
            if self.config.verbose:
                self.log.info(line)
            results.add_seg_result(amp, 'READ_NOISE', Nread)
        results.write(clobber=True)
