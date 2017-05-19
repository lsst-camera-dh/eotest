"""
@brief Compute detector response vs incident flux for pairs of flats.
These data are to be used for linearity and full-well measurments.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import glob
import numpy as np
import astropy.io.fits as fits
from lsst.eotest.fitsTools import fitsTableFactory, fitsWriteto
import lsst.eotest.image_utils as imutils
from .MaskedCCD import MaskedCCD
from .EOTestResults import EOTestResults
from .DetectorResponse import DetectorResponse

import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

def pair_mean(flat1, flat2, amp):
    """
    Compute the mean pixel value for a given segment, averaged between
    two exposures (assumed to have the same exposure time). The inputs
    flat1 and flat2 are MaskedCCD objects.  Any masks are assumed to
    have been enabled in the stat_ctrl attribute.
    """
    #
    # Remove bias using fit to overscan region and trim to imaging region.
    #
    im1 = flat1.unbiased_and_trimmed_image(amp)
    im2 = flat2.unbiased_and_trimmed_image(amp)
    #
    # Compute the mean DN of the images.
    #
    stats1 = afwMath.makeStatistics(im1, afwMath.MEAN, flat1.stat_ctrl)
    stats2 = afwMath.makeStatistics(im2, afwMath.MEAN, flat2.stat_ctrl)

    avg_mean_value = (stats1.getValue(afwMath.MEAN) +
                      stats2.getValue(afwMath.MEAN))/2.
    return avg_mean_value

def find_flat2(flat1):
    pattern = flat1.split('flat1')[0] + 'flat2*.fits'
    flat2 = glob.glob(pattern)[0]
    return flat2

class FlatPairConfig(pexConfig.Config):
    """Configuration for flat pair task"""
    output_dir = pexConfig.Field("Output directory", str, default='.')
    eotest_results_file = pexConfig.Field("EO test results filename",
                                          str, default=None)
    verbose = pexConfig.Field("Turn verbosity on", bool, default=True)

class FlatPairTask(pipeBase.Task):
    """Task to compute detector response vs incident flux from
       flat pair dataset."""
    ConfigClass = FlatPairConfig
    _DefaultName = "FlatPairTask"

    @pipeBase.timeMethod
    def run(self, sensor_id, infiles, mask_files, gains, detrespfile=None,
            bias_frame=None, max_pd_frac_dev=0.05,
            linearity_spec_range=(1e3, 9e4), use_exptime=False):
        self.sensor_id = sensor_id
        self.infiles = infiles
        self.mask_files = mask_files
        self.gains = gains
        self.bias_frame = bias_frame
        self.max_pd_frac_dev = max_pd_frac_dev
        if detrespfile is None:
            #
            # Compute detector response from flat pair files.
            #
            detrespfile = self.extract_det_response(use_exptime)
        #
        # Perform full well and linearity analyses.
        #
        detresp = DetectorResponse(detrespfile)

        outfile = self.config.eotest_results_file
        if outfile is None:
            outfile = os.path.join(self.config.output_dir,
                                   '%s_eotest_results.fits' % self.sensor_id)
        all_amps = imutils.allAmps(detrespfile)
        output = EOTestResults(outfile, namps=len(all_amps))
        if self.config.verbose:
            self.log.info("Amp        full well (e-/pixel)   max. frac. dev.")
        for amp in all_amps:
            try:
                full_well, fp = detresp.full_well(amp)
            except StandardError as eobj:
                self.log.info("Exception caught in full well calculation:")
                self.log.info(str(eobj))
                full_well = None
            self.log.info('linearity analysis range: %s, %s' %
                          linearity_spec_range)
            try:
                maxdev, fit_pars, Ne, flux = \
                    detresp.linearity(amp, spec_range=linearity_spec_range)
            except StandardError as eobj:
                self.log.info("Exception caught in linearity calculation:")
                self.log.info(str(eobj))
                maxdev = None
            if self.config.verbose:
                self.log.info('%2i            %s             %s'
                              % (amp, full_well, maxdev))
            if full_well is not None:
                output.add_seg_result(amp, 'FULL_WELL', full_well)
            if maxdev is not None:
                output.add_seg_result(amp, 'MAX_FRAC_DEV', float(maxdev))
        output.write()
    def _create_detresp_fits_output(self, nrows):
        self.output = fits.HDUList()
        self.output.append(fits.PrimaryHDU())
        all_amps = imutils.allAmps()
        colnames = ['flux'] + ['AMP%02i_SIGNAL' % i for i in all_amps]
        formats = 'E'*len(colnames)
        units = ['None'] + ['e-']*len(all_amps)
        columns = [np.zeros(nrows, dtype=np.float) for fmt in formats]
        fits_cols = [fits.Column(name=colnames[i], format=formats[i],
                                 unit=units[i], array=columns[i])
                     for i in range(len(units))]
        hdu = fitsTableFactory(fits_cols)
        hdu.name = 'DETECTOR_RESPONSE'
        self.output.append(hdu)
    def extract_det_response(self, use_exptime):
        max_pd_frac_dev = self.max_pd_frac_dev
        outfile = os.path.join(self.config.output_dir,
                               '%s_det_response.fits' % self.sensor_id)
        file1s = sorted([item for item in self.infiles
                         if item.find('flat1') != -1])
        if self.config.verbose:
            self.log.info("writing to %s" % outfile)
        self._create_detresp_fits_output(len(file1s))
        for row, file1 in enumerate(file1s):
            try:
                file2 = find_flat2(file1)
            except IndexError:
                # Just use flat1 again since only average is taken and
                # FPN subtraction isn't needed.
                file2 = file1

            if self.config.verbose:
                self.log.info("processing\n   %s as flat1 and\n   %s as flat2"
                              % (file1, file2))

            flat1 = MaskedCCD(file1, mask_files=self.mask_files,
                              bias_frame=self.bias_frame)
            flat2 = MaskedCCD(file2, mask_files=self.mask_files,
                              bias_frame=self.bias_frame)

            pd1 = flat1.md.get('MONDIODE')
            pd2 = flat2.md.get('MONDIODE')
            exptime1 = flat1.md.get('EXPTIME')
            exptime2 = flat2.md.get('EXPTIME')

            if exptime1 != exptime2:
                raise RuntimeError("Exposure times do not match for:\n%s\n%s\n"
                                   % (file1, file2))
            if (not use_exptime or
                ((type(pd1) != str and type(pd2) != str) and
                 (pd1 != 0 and pd2 != 0))):
                flux = abs(pd1*exptime1 + pd2*exptime2)/2.
                if np.abs((pd1 - pd2)/((pd1 + pd2)/2.)) > max_pd_frac_dev:
                    self.log.info("Skipping %s and %s since MONDIODE values do not agree to %.1f%%" % (file1, file2, max_pd_frac_dev*100.))
                    continue
            else:
                flux = exptime1
            if self.config.verbose:
                self.log.info('   row = %s' % row)
                self.log.info('   pd1, pd2 = %s, %s' % (pd1, pd2))
                self.log.info('   exptime1, exptime2 = %s, %s '
                              % (exptime1, exptime2))
                self.log.info('   flux = %s' % flux)
                self.log.info('   flux/exptime = %s' % (flux/exptime1,))
            self.output[-1].data.field('FLUX')[row] = flux
            for amp in flat1:
                # Convert to e- and write out for each segment.
                signal = pair_mean(flat1, flat2, amp)*self.gains[amp]
                self.output[-1].data.field('AMP%02i_SIGNAL' % amp)[row] = signal
        self.output[0].header['NAMPS'] = len(flat1)
        fitsWriteto(self.output, outfile, clobber=True)
        return outfile
