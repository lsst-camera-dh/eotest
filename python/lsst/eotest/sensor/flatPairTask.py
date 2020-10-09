"""
@brief Compute detector response vs incident flux for pairs of flats.
These data are to be used for linearity and full-well measurments.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import glob
import numpy as np
import astropy.io.fits as fits
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.eotest.fitsTools import fitsTableFactory, fitsWriteto
import lsst.eotest.image_utils as imutils
from .MaskedCCD import MaskedCCD
from .EOTestResults import EOTestResults
from .DetectorResponse import DetectorResponse


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

    flat1_value = stats1.getValue(afwMath.MEAN)
    flat2_value = stats2.getValue(afwMath.MEAN)
    avg_mean_value = (flat1_value + flat2_value)/2.
    return np.array([avg_mean_value, flat1_value, flat2_value], dtype=float)


def row_mean_variance(flat1, flat2, amp):
    """
    Compute the 3-sigma clipped variance of the row mean distributions
    for the specified amp from a pair of flats.
    """
    mi_diff = afwImage.MaskedImageF(flat1.unbiased_and_trimmed_image(amp),
                                    deep=True)
    mi_diff -= flat2.unbiased_and_trimmed_image(amp)
    row_means = np.mean(mi_diff.getImage().array, axis=1)
    return afwMath.makeStatistics(row_means, afwMath.VARIANCECLIP).getValue()


def find_flat2(flat1):
    pattern = flat1.split('flat1')[0] + 'flat2*.fits'
    flat2_files = glob.glob(pattern)
    if len(flat2_files) == 1:
        return flat2_files[0]
    with fits.open(flat1) as hdus:
        exptime1 = hdus[0].header['EXPTIME']
    for flat2 in flat2_files:
        with fits.open(flat2) as hdus:
            if hdus[0].header['EXPTIME'] == exptime1:
                return flat2
    raise RuntimeError("no flat2 file found for {}".format(flat1))


def mondiode_value(fits_file, _, factor=5):
    """Compute the effective monitoring diode current by integrating
    over the pd current time history in the AMP0_MEAS_TIMES extension
    and dividing by the EXPTIME value.
    """
    with fits.open(fits_file) as hdus:
        x = hdus['AMP0.MEAS_TIMES'].data.field('AMP0_MEAS_TIMES')
        y = -hdus['AMP0.MEAS_TIMES'].data.field('AMP0_A_CURRENT')*1e9
        exptime = hdus[0].header['EXPTIME']
    ythresh = (max(y) - min(y))/factor + min(y)
    index = np.where(y < ythresh)
    y0 = np.median(y[index])
    y -= y0
    return sum((y[1:] + y[:-1])/2.*(x[1:] - x[:-1]))/exptime


def compute_row_mean_var_slopes(detrespfile, min_flux=3000, max_flux=1e5):
    """
    Fits linear slopes to var(row_mean) vs flux as an indicator of
    long-range serial correlations.
    """
    slopes = dict()
    with fits.open(detrespfile) as detresp:
        ncols = detresp[0].header['NUMCOLS']
        amps = range(1, detresp[0].header['NAMPS'] + 1)
        for amp in amps:
            amp_label = f'AMP{amp:02d}'
            flux = detresp[1].data[f'{amp_label}_SIGNAL']
            row_mean_var = detresp[1].data[f'{amp_label}_ROW_MEAN_VAR']
            # Restrict to higher flux values below full well and
            # avoid nans in row_mean_var.
            index = np.where((min_flux < flux) & (flux < max_flux)
                             & (row_mean_var == row_mean_var))
            if len(index[0]) == 0:
                slopes[amp] = 0
            else:
                slopes[amp] = sum(row_mean_var[index])/sum(2.*flux[index]/ncols)
    return slopes


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
            linearity_spec_range=(1e3, 9e4), use_exptime=False,
            flat2_finder=find_flat2, mondiode_func=mondiode_value,
            linearity_correction=None, dark_frame=None):
        self.sensor_id = sensor_id
        self.infiles = infiles
        self.mask_files = mask_files
        self.gains = gains
        self.bias_frame = bias_frame
        self.max_pd_frac_dev = max_pd_frac_dev
        self.find_flat2 = flat2_finder
        self.mondiode_func = mondiode_func
        self.linearity_correction = linearity_correction
        self.dark_frame = dark_frame
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
        row_mean_var_slopes = compute_row_mean_var_slopes(detrespfile)
        if self.config.verbose:
            self.log.info("Amp       max_signal (e-/pixel)   max. frac. dev.")
        for amp in all_amps:
            try:
                full_well, _ = detresp.full_well(amp)
            except Exception as eobj:
                self.log.info("Exception caught in full well calculation:")
                self.log.info(str(eobj))
                full_well = None
            self.log.info('linearity analysis range: %s, %s' %
                          linearity_spec_range)
            try:
                results = \
                    detresp.linearity(amp, spec_range=linearity_spec_range)
                maxdev = results[0]
                turnoff = results[-1]
            except Exception as eobj:
                self.log.info("Exception caught in linearity calculation:")
                self.log.info(str(eobj))
                maxdev = None
                turnoff = None
            # The maximum observed signal should be reported in ADUs.
            max_observed_signal = np.max(detresp.Ne[amp])/self.gains[amp]
            if self.config.verbose:
                self.log.info('%2i            %s             %s'
                              % (amp, max_observed_signal, maxdev))
            if full_well is not None:
                output.add_seg_result(amp, 'FULL_WELL', full_well)
            if maxdev is not None:
                output.add_seg_result(amp, 'MAX_FRAC_DEV', float(maxdev))
            output.add_seg_result(amp, 'ROW_MEAN_VAR_SLOPE',
                                  row_mean_var_slopes[amp])
            output.add_seg_result(amp, 'MAX_OBSERVED_SIGNAL',
                                  max_observed_signal)
            if turnoff is not None:
                output.add_seg_result(amp, 'LINEARITY_TURNOFF', float(turnoff))
        output.write()

    def _create_detresp_fits_output(self, nrows, infile):
        self.output = fits.HDUList()
        self.output.append(fits.PrimaryHDU())
        all_amps = imutils.allAmps(infile)
        colnames = ['flux'] + ['AMP%02i_SIGNAL' % i for i in all_amps] + \
                   ['FLAT1_AMP%02i_SIGNAL' % i for i in all_amps] + \
                   ['FLAT2_AMP%02i_SIGNAL' % i for i in all_amps] + \
                   ['AMP%02i_ROW_MEAN_VAR' % i for i in all_amps] + \
                   ['SEQNUM', 'DAYOBS']
        formats = 'E'*(len(colnames) - 2) + 'JJ'
        namps = len(all_amps)
        units = ['None'] + 3*namps*['e-'] + namps*['e-^2'] + ['None', 'None']
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
        self._create_detresp_fits_output(len(file1s), file1s[0])
        for row, file1 in enumerate(file1s):
            try:
                file2 = self.find_flat2(file1)
            except IndexError:
                # Just use flat1 again since only average is taken and
                # FPN subtraction isn't needed.
                file2 = file1

            if self.config.verbose:
                self.log.info("processing\n   %s as flat1 and\n   %s as flat2",
                              file1, file2)

            flat1 = MaskedCCD(file1, mask_files=self.mask_files,
                              bias_frame=self.bias_frame,
                              dark_frame=self.dark_frame,
                              linearity_correction=self.linearity_correction)
            flat2 = MaskedCCD(file2, mask_files=self.mask_files,
                              bias_frame=self.bias_frame,
                              dark_frame=self.dark_frame,
                              linearity_correction=self.linearity_correction)

            seqnum = flat1.md.get('SEQNUM')
            try:
                dayobs = flat1.md.get('DAYOBS')
            except KeyError:
                # This occurs if running on non-BOT data.
                dayobs = 0

            exptime1 = flat1.md.get('EXPTIME')
            exptime2 = flat2.md.get('EXPTIME')

            if exptime1 != exptime2:
                raise RuntimeError("Exposure times do not match for:\n%s\n%s\n"
                                   % (file1, file2))

            if self.mondiode_func is None:
                pd1 = flat1.md.get('MONDIODE')
                pd2 = flat2.md.get('MONDIODE')
            else:
                try:
                    pd1 = self.mondiode_func(file1, exptime1)
                    pd2 = self.mondiode_func(file2, exptime2)
                except KeyError as eobj:
                    self.log.info("KeyError exception computing pd current:\n"
                                  + str(eobj))
                    continue

            if use_exptime:
                flux = exptime1
            else:
                flux = abs(pd1*exptime1 + pd2*exptime2)/2.
                if np.abs((pd1 - pd2)/((pd1 + pd2)/2.)) > max_pd_frac_dev:
                    self.log.info("Skipping %s and %s since MONDIODE values "
                                  "do not agree to %.1f%%",
                                  file1, file2, max_pd_frac_dev*100.)
                    continue
            if self.config.verbose:
                self.log.info('   row = %s', row)
                self.log.info('   pd1, pd2 = %s, %s', pd1, pd2)
                self.log.info('   exptime1, exptime2 = %s, %s ',
                              exptime1, exptime2)
                self.log.info('   flux = %s', flux)
                self.log.info('   flux/exptime = %s', flux/exptime1)
            self.output[-1].data.field('FLUX')[row] = flux
            for amp in flat1:
                # Convert to e- and write out for each segment.
                signal, sig1, sig2 \
                    = pair_mean(flat1, flat2, amp)*self.gains[amp]
                colname = 'AMP%02i_SIGNAL' % amp
                self.output[-1].data.field(colname)[row] = signal
                self.output[-1].data.field('FLAT1_' + colname)[row] = sig1
                self.output[-1].data.field('FLAT2_' + colname)[row] = sig2
                self.output[-1].data.field(f'AMP{amp:02d}_ROW_MEAN_VAR')[row] \
                    = row_mean_variance(flat1, flat2, amp)
                self.output[-1].data.field('SEQNUM')[row] = seqnum
                self.output[-1].data.field('DAYOBS')[row] = dayobs
        self.output[0].header['NAMPS'] = len(flat1)
        self.output[0].header['NUMCOLS'] = flat1.amp_geom.imaging.getWidth()
        fitsWriteto(self.output, outfile, overwrite=True)
        return outfile
