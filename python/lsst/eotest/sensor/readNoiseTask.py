"""
@brief Compute read noise distributions for a sample of images.  Bias
and system readout noise exposures, the latter for determining the
noise contribution from the electronics, must be provided.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import numpy as np
import astropy.io.fits as fits
from lsst.eotest.fitsTools import fitsTableFactory, fitsWriteto
import lsst.eotest.image_utils as imutils
from MaskedCCD import MaskedCCD
from EOTestResults import EOTestResults
from read_noise import noise_dists, NoiseDistributions
import lsst.afw.geom as afwGeom

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

def _write_read_noise_dists(outfile, Ntot, Nsys, gains, bias, sysnoise):
    output = fits.HDUList()
    output.append(fits.PrimaryHDU())
    output[0].header['BIASFILE'] = bias
    output[0].header['SYSNFILE'] = str(sysnoise)
    for amp in Ntot:
        sigtot, sigsys = Ntot[amp], Nsys[amp]
        nread_col = fits.Column(name="TOTAL_NOISE", format="E",
                                unit="e- rms", array=sigtot)
        nsys_col = fits.Column(name="SYSTEM_NOISE", format="E",
                               unit="e- rms", array=sigsys)
        output.append(fitsTableFactory((nread_col, nsys_col)))
        output[amp].name = "SEGMENT%s" % imutils.channelIds[amp]
        output[0].header["GAIN%s" % imutils.channelIds[amp]] = gains[amp]
        output[0].header["SIGTOT%s" % imutils.channelIds[amp]] = \
            imutils.median(sigtot)
        output[0].header["SIGSYS%s" % imutils.channelIds[amp]] = \
            imutils.median(sigsys)
    fitsWriteto(output, outfile, clobber=True)

class ReadNoiseConfig(pexConfig.Config):
    """Configuration for read noise task"""
    dx = pexConfig.Field("Subregion size in pixels along x-direction",
                         int, default=100)
    dy = pexConfig.Field("Subregion size in pixels along y-direction",
                         int, default=100)
    nsamp = pexConfig.Field("Number of subregions to sample",
                            int, default=1000)
    temp_set_point = pexConfig.Field("Required temperature (C) set point",
                                     float, default=-95.)
    temp_set_point_tol = pexConfig.Field("Required temperature set point tolerance (degrees C)",
                                         float, default=1.)
    output_dir = pexConfig.Field("Output directory", str, default=".")
    eotest_results_file = pexConfig.Field("EO test results filename",
                                          str, default=None)
    verbose = pexConfig.Field("Turn verbosity on", bool, default=True)

class ReadNoiseTask(pipeBase.Task):
    """Task to estimate sensor read noise."""
    ConfigClass = ReadNoiseConfig
    _DefaultName = "ReadNoiseTask"

    @pipeBase.timeMethod
    def run(self, sensor_id, bias_files, gains, system_noise_files=None,
            system_noise=None, mask_files=(), use_overscan=False):
        imutils.check_temperatures(bias_files, self.config.temp_set_point_tol,
                                   setpoint=self.config.temp_set_point,
                                   warn_only=True)
        outfiles = []
        Ntot = NoiseDistributions()
        Nsys = NoiseDistributions()
        if system_noise_files is None:
            system_noise_files = [None]*len(bias_files)
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
            if use_overscan:
                imaging = ccd.amp_geom.serial_overscan
                dx = imaging.getWidth()/2
                dy = self.config.dy
                nsamp = self.config.nsamp
            else:
                imaging = ccd.amp_geom.imaging
                dx = self.config.dx
                dy = self.config.dy
                nsamp = self.config.nsamp
            #
            # Create a single sub-region sampler so that the same
            # sub-regions will be used for both the bias and system
            # noise frames.
            #
            sampler = imutils.SubRegionSampler(dx, dy, nsamp, imaging=imaging)

            Ntot_amp = noise_dists(bias, gains, sampler, mask_files=mask_files)
            Nsys_amp = noise_dists(sysnoise, gains, sampler,
                                   mask_files=mask_files)

            _write_read_noise_dists(outfile, Ntot_amp, Nsys_amp, gains,
                                    bias, sysnoise)
            #
            # Accumulate noise distributions for final median calculation
            #
            Ntot.append(Ntot_amp)
            Nsys.append(Nsys_amp)

        results_file = self.config.eotest_results_file
        if results_file is None:
            results_file = os.path.join(self.config.output_dir,
                                        '%s_eotest_results.fits' % sensor_id)

        results = EOTestResults(results_file)
        if self.config.verbose:
            self.log.info("Amp    read noise    total noise    system noise")
        for amp in ccd:
            Ntot_med = imutils.median(Ntot[amp])
            if system_noise is not None:
                Nsys_med = float(system_noise[amp])
            else:
                Nsys_med = imutils.median(Nsys[amp])
            var = Ntot_med**2 - Nsys_med**2
            if var >= 0:
                Nread = np.sqrt(var)
            else:
                Nread = -1
            line = "%2s       %7.4f        %7.4f        %7.4f" % (amp, Nread,
                                                                  Ntot_med,
                                                                  Nsys_med)
            if self.config.verbose:
                self.log.info(line)
            results.add_seg_result(amp, 'READ_NOISE', Nread)
            results.add_seg_result(amp, 'TOTAL_NOISE', Ntot_med)
            results.add_seg_result(amp, 'SYSTEM_NOISE', Nsys_med)
        results.write(clobber=True)
