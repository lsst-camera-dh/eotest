"""
@brief For pairs of flats obtain for a range of exposures, compute the
photon transfer curve and compute and write out the full well.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
from __future__ import print_function
import os
import glob
import operator
import numpy as np
import scipy.optimize
import astropy.io.fits as fits
from lsst.eotest.fitsTools import fitsTableFactory, fitsWriteto
import lsst.eotest.image_utils as imutils
from .MaskedCCD import MaskedCCD
from .EOTestResults import EOTestResults

import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase


def find_flat2(flat1):
    pattern = flat1.split('flat1')[0] + 'flat2*.fits'
    try:
        flat2 = glob.glob(pattern)[0]
        return flat2
    except IndexError:
        return flat1


def exptime(x): return imutils.Metadata(x).get('EXPTIME')


def glob_flats(full_path, outfile='ptc_flats.txt'):
    flats = glob.glob(os.path.join(full_path, '*_flat?.fits'))
    output = open(outfile, 'w')
    for item in flats:
        output.write('%s\n' % item)
    output.close()


def find_flats(args):
    files = args.files(args.flats, args.flats_file_list)
    file1s = sorted([item.strip() for item in files
                     if item.find('flat1') != -1])
    return [(f1, find_flat2(f1)) for f1 in file1s]


def ptc_func(pars, mean):
    """
    Model for variance vs mean.
    See http://adsabs.harvard.edu/abs/2015A%26A...575A..41G.
    """
    alpha, gain = pars
    return mean*(1./gain - mean*alpha)


def residuals(pars, mean, var):
    "Residuals function for least-squares fit of PTC curve."
    return (var - ptc_func(pars, mean))/np.sqrt(var)


class FlatPairStats(object):
    def __init__(self, fmean, fvar):
        self.flat_mean = fmean
        self.flat_var = fvar

def flat_pair_stats(ccd1, ccd2, amp, mask_files=(), bias_frame=None):
    if ccd1.md.get('EXPTIME') != ccd2.md.get('EXPTIME'):
        raise RuntimeError("Exposure times for files %s, %s do not match"
                           % (ccd1.imfile, ccd2.imfile))
    #
    # Mean and variance calculations that account for masks (via
    # ccd1.stat_ctrl, which is the same for both MaskedImages).
    #

    def mean(im): return afwMath.makeStatistics(im, afwMath.MEAN,
                                                ccd1.stat_ctrl).getValue()

    def var(im): return afwMath.makeStatistics(im, afwMath.VARIANCE,
                                               ccd1.stat_ctrl).getValue()
    #
    # Extract imaging region for segments of both CCDs.
    #
    image1 = ccd1.unbiased_and_trimmed_image(amp, bias_frame=bias_frame)
    image2 = ccd2.unbiased_and_trimmed_image(amp, bias_frame=bias_frame)
    if ccd1.imfile == ccd2.imfile:
        # Don't have pairs of flats, so estimate noise and gain
        # from a single frame, ignoring FPN.
        fmean = mean(image1)
        fvar = var(image1)
    else:
        #
        # Make a deep copy since otherwise the pixel values in image1
        # would be altered in the ratio calculation.
        #
        fratio_im = afwImage.MaskedImageF(image1, True)
        operator.itruediv(fratio_im, image2)
        fratio = mean(fratio_im)
        image2 *= fratio
        fmean = (mean(image1) + mean(image2))/2.

        fdiff = afwImage.MaskedImageF(image1, True)
        fdiff -= image2
        fvar = var(fdiff)/2.

    return FlatPairStats(fmean, fvar)


class PtcConfig(pexConfig.Config):
    """Configuration for ptc task"""
    output_dir = pexConfig.Field("Output directory", str, default='.')
    eotest_results_file = pexConfig.Field("EO test results filename",
                                          str, default=None)
    max_frac_offset = pexConfig.Field(
        "maximum fraction offset from median gain curve to omit points from PTC fit.", float, default=0.2)
    verbose = pexConfig.Field("Turn verbosity on", bool, default=True)


class PtcTask(pipeBase.Task):
    """Task to compute photon transfer curve from flat pair dataset"""
    ConfigClass = PtcConfig
    _DefaultName = "PtcTask"

    @pipeBase.timeMethod
    def run(self, sensor_id, infiles, mask_files, gains, binsize=1,
            bias_frame=None):
        outfile = os.path.join(self.config.output_dir,
                               '%s_ptc.fits' % sensor_id)
        all_amps = imutils.allAmps(infiles[0])
        ptc_stats = dict([(amp, ([], [])) for amp in all_amps])
        exposure = []
        file1s = sorted([item for item in infiles if item.find('flat1') != -1])
        for flat1 in file1s:
            flat2 = find_flat2(flat1)
            if self.config.verbose:
                self.log.info("processing %s" % flat1)
            exposure.append(exptime(flat1))
            ccd1 = MaskedCCD(flat1, mask_files=mask_files,
                             bias_frame=bias_frame)
            ccd2 = MaskedCCD(flat2, mask_files=mask_files,
                             bias_frame=bias_frame)
            for amp in ccd1:
                results = flat_pair_stats(ccd1, ccd2, amp,
                                          mask_files=mask_files,
                                          bias_frame=bias_frame) 
                ptc_stats[amp][0].append(results.flat_mean)
                ptc_stats[amp][1].append(results.flat_var)
        self._fit_curves(ptc_stats, sensor_id)
        output = fits.HDUList()
        output.append(fits.PrimaryHDU())
        colnames = ['EXPOSURE']
        units = ['seconds']
        columns = [np.array(exposure, dtype=np.float)]
        for amp in all_amps:
            colnames.extend(['AMP%02i_MEAN' % amp, 'AMP%02i_VAR' % amp])
            units.extend(['ADU', 'ADU**2'])
            columns.extend([np.array(ptc_stats[amp][0], dtype=np.float),
                            np.array(ptc_stats[amp][1], dtype=np.float)])
        formats = 'E'*len(colnames)
        fits_cols = [fits.Column(name=colnames[i], format=formats[i],
                                 unit=units[i], array=columns[i])
                     for i in range(len(columns))]
        output.append(fitsTableFactory(fits_cols))
        output[-1].name = 'PTC_STATS'
        output[0].header['NAMPS'] = len(all_amps)
        fitsWriteto(output, outfile, overwrite=True)

    def _fit_curves(self, ptc_stats, sensor_id):
        """
        Fit a quadratic to the variance vs mean data for each amp.
        See http://adsabs.harvard.edu/abs/2015A%26A...575A..41G.
        """
        outfile = self.config.eotest_results_file
        if outfile is None:
            outfile = os.path.join(self.config.output_dir,
                                   '%s_eotest_results.fits' % sensor_id)
        output = EOTestResults(outfile, namps=len(ptc_stats))
        for amp in ptc_stats:
            mean, var = np.array(ptc_stats[amp][0]), np.array(ptc_stats[amp][1])
            # Compute the median gain and remove outliers.
            med_gain = np.median(mean/var)
            frac_resids = np.abs((var - mean/med_gain)/var)
            index = np.where(frac_resids < self.config.max_frac_offset)
            # Least-squares fit to the brighter-fatter model.
            p0 = 0, med_gain
            try:
                results = scipy.optimize.leastsq(residuals, p0, full_output=1,
                                                 args=(mean[index], var[index]))
                pars, cov = results[:2]
                ptc_gain = pars[1]
                ptc_error = np.sqrt(cov[1][1])
            except Exception as eobj:
                self.log.info("Exception caught while fitting PTC:")
                self.log.info(str(eobj))
                ptc_gain = 0.
                ptc_error = -1.
            # Write gain and error to EO test results file.
            output.add_seg_result(amp, 'PTC_GAIN', ptc_gain)
            output.add_seg_result(amp, 'PTC_GAIN_ERROR', ptc_error)
            self.log.info("%i  %f  %f" % (amp, ptc_gain, ptc_error))
        output.write()
