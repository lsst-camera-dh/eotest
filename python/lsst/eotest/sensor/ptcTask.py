"""
@brief For pairs of flats obtain for a range of exposures, compute the
photon transfer curve and compute and write out the full well.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import glob
from copy import deepcopy
import operator
import numpy as np
import scipy.optimize
import astropy.io.fits as fits
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.eotest.fitsTools import fitsTableFactory, fitsWriteto
import lsst.eotest.image_utils as imutils
from .MaskedCCD import MaskedCCD
from .EOTestResults import EOTestResults
from .flatPairTask import find_flat2
import astropy.stats as astats
import matplotlib.pyplot as plt

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
    Model for variance vs mean.  See Astier et al.
    https://confluence.slac.stanford.edu/pages/viewpage.action?pageId=242286867
    """
    a00, gain, intcpt = pars
    return 0.5/(a00*gain*gain)*(1 - np.exp(-2*a00*mean*gain)) + intcpt/(gain*gain)


def residuals(pars, mean, var):
    """
    Residuals function for least-squares fit of PTC curve.
    """
    return (var - ptc_func(pars, mean))/np.sqrt(var)


class FlatPairStats(object):
    def __init__(self, fmean, fvar, fmean_mad, fvar_mad):
        self.flat_mean = fmean
        self.flat_var = fvar
        self.flat_mean_mad = fmean_mad
        self.flat_var_mad = fvar_mad

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
        #fratio_im = afwImage.MaskedImageF(image1, True)
        #operator.itruediv(fratio_im, image2)
        #fratio1 = mean(fratio_im)
        fratio1 = mean(image1)/mean(image2)
        #print('Before scaling')
        #print('mean(image2):  ' + str(mean(image2)))
        image2 *= fratio1
        #print('After scaling')
        #print('mean(image2):  ' + str(mean(image2)))
        fmean = (mean(image1) + mean(image2))/2.

        fdiff1 = afwImage.MaskedImageF(image1, True)
        fdiff1 -= image2
        #print(amp,mean(image1),mean(image2),fmean, mean(fdiff))
        fvar = var(fdiff1)/2.

        #print(np.min(np.ravel(fdiff.getArrays()[0])), np.max(np.ravel(fdiff.getArrays()[0])))
        #if amp == 7:
            #fig,ax = plt.subplots()
            #temp = plt.hist(np.ravel(fdiff.getArrays()[0]),bins=200,histtype='step',color='orange')
            #plt.plot([-np.sqrt(fvar), -np.sqrt(fvar)],[0,np.max(temp[0])],color='blue')
            #plt.plot([np.sqrt(fvar), np.sqrt(fvar)],[0,np.max(temp[0])],color='blue')
        
	# Make a more robust estimate of variance
        fdiff = np.ravel(fdiff1.getArrays()[0])
        fvar_mad1 = astats.mad_std(fdiff)  #/2.
        image1 = np.ravel(image1.getArrays()[0])
        #image2 = np.ravel(ccd2.unbiased_and_trimmed_image(amp, bias_frame=bias_frame).getArrays()[0])
        image2 = np.ravel(image2.getArrays()[0])
        #g = np.where((np.abs(fdiff) < (fvar_mad1*14.826))*(np.abs(image2) > 1e-3))[0]
        #g1 = np.where((np.abs(image2) > 1e-3))[0]
        g = np.where((np.abs(fdiff) < (fvar_mad1*14.826)))[0]
        #if amp == 15:
            #fig,ax = plt.subplots()
            #plt.hist(fdiff,bins=200)
            #print('Mean image1:  ' + str(mean(image1)))
            #print('Mean image2:  ' + str(mean(image2)))

        #fratio = np.mean(image1[g]/image2[g])
        fratio = np.mean(image1[g])/np.mean(image2[g])
        image2 *= fratio
        fmean_mad = (np.mean(image1[g]) + np.mean(image2[g]))/2.
        fdiff = image1[g] - image2[g]
        fvar_mad = np.var(fdiff)/2.
        print('len(g):  ' + str(len(g)))
        if fmean_mad < 0:
            print('!!!!!!!!!!! fmean_mad:  ' + str(fmean_mad))
            print('len(g):  ' + str(len(g)))
            #print('len(g1):  ' + str(len(g1)))
            #print('len(g2):  ' + str(len(g2)))
            print('fmean:  ' + str(fmean))
            print('min,max fdiff:  ' + str(min(fdiff)) + ' ' + str(max(fdiff)))
            print('min,max fdiff1:  ' + str(np.min(fdiff1.getArrays()[0])) + ' ' + str(np.max(fdiff1.getArrays()[0])))
            print('str(type(fdiff1)):  ' + str(type(fdiff1)))
            print('fmean_mad:  ' + str(fmean_mad))
            print('np.mean(image1[g]):  ' + str(np.mean(image1[g])))
            print('np.mean(image2[g]):  ' + str(np.mean(image2[g])))
            print('np.mean(image1):  ' + str(np.mean(image1)))
            print('np.mean(image2):  ' + str(np.mean(image2)))
            print('fratio:  ' + str(fratio))
            print('fratio1:  ' + str(fratio1))
            print('fvar:  ' + str(fvar))
            print('fvar_mad:  ' + str(fvar_mad))
            print('fvar_mad1:  ' + str(fvar_mad1))
            print('amp:  ' + str(amp))
            print('mean(fdiff1):  ' + str(mean(fdiff1)))
            print('mean(fdiff):  ' + str(mean(fdiff)))
            print('files:  ')
            print(ccd1.imfile)
            print(ccd2.imfile)


        #if amp == 7:
            #temp = plt.hist(np.ravel(fdiff),bins=200,color='green',histtype='step')
            #plt.plot([-np.sqrt(fvar_mad), -np.sqrt(fvar_mad)],[0,np.max(temp[0])],color='red')
            #plt.plot([np.sqrt(fvar_mad), np.sqrt(fvar_mad)],[0,np.max(temp[0])],color='red')
            #plt.title('AMP' + str(amp).zfill(2))
            #outfile = 'temp_png/hist_' + ccd1.imfile.split('/')[-1].split('.')[0] + '.png'
            #print('Saving image file:  ' + outfile)
            #plt.savefig(outfile)

        #print('fmean:  ' + str(fmean))
        #print('fvar:  ' + str(fvar))
        #print('fmean_mad:  ' + str(fmean_mad))
        #print('fvar_mad:  ' + str(fvar_mad))
    return FlatPairStats(fmean, fvar, fmean_mad, fvar_mad)


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
            bias_frame=None, flat2_finder=find_flat2,
            linearity_correction=None):
        outfile = os.path.join(self.config.output_dir,
                               '%s_ptc.fits' % sensor_id)
        all_amps = imutils.allAmps(infiles[0])
        #print(all_amps)
        ptc_stats = dict([(amp, ([], [], [], [])) for amp in all_amps])
        exposure = []
        seqnums = []
        dayobs = []
        file1s = sorted([item for item in infiles if item.find('flat0') != -1])
        print('# file1s:  ' + str(len(file1s)))
        for flat1 in file1s:
            flat2 = flat2_finder(flat1)
            if self.config.verbose:
                self.log.info("processing %s" % flat1)
            exposure.append(exptime(flat1))
            print('Reading ccd1')
            print(flat1)
            ccd1 = MaskedCCD(flat1, mask_files=mask_files,
                             bias_frame=bias_frame,
                             linearity_correction=linearity_correction)
            #ccd1.stat_ctrl.setNumIter(3)
            #ccd1.stat_ctrl.setNumSigmaClip(10.0)
            print('Reading ccd2')
            print(flat2)
            ccd2 = MaskedCCD(flat2, mask_files=mask_files,
                             bias_frame=bias_frame,
                             linearity_correction=linearity_correction)
            #ccd2.stat_ctrl.setNumIter(3)
            #ccd2.stat_ctrl.setNumSigmaClip(10.0)
            #print(len(ccd1), len(ccd2))
            for amp in ccd1:
                results = flat_pair_stats(ccd1, ccd2, amp,
                                          mask_files=mask_files,
                                          bias_frame=bias_frame)
                ptc_stats[amp][0].append(results.flat_mean)
                ptc_stats[amp][1].append(results.flat_var)
                ptc_stats[amp][2].append(results.flat_mean_mad)
                ptc_stats[amp][3].append(results.flat_var_mad)
            seqnums.append(ccd1.md.get('SEQNUM'))
            try:
                dayobs.append(ccd1.md.get('DAYOBS'))
            except KeyError:
                dayobs.append(0)
        print('Now fitting the curves')
        self._fit_curves(ptc_stats, sensor_id)
        output = fits.HDUList()
        output.append(fits.PrimaryHDU())
        colnames = ['EXPOSURE']
        units = ['seconds']
        columns = [np.array(exposure, dtype=np.float)]
        for amp in all_amps:
            colnames.extend(['AMP%02i_MEAN' % amp, 'AMP%02i_VAR' % amp, 'AMP%02i_MMEAN' % amp, 'AMP%02i_MVAR' % amp])
            units.extend(['ADU', 'ADU**2', 'ADU', 'ADU**2'])
            columns.extend([np.array(ptc_stats[amp][0], dtype=np.float),
                            np.array(ptc_stats[amp][1], dtype=np.float),
                            np.array(ptc_stats[amp][2], dtype=np.float),
                            np.array(ptc_stats[amp][3], dtype=np.float)])
        colnames.append('SEQNUM')
        units.append('None')
        columns.append(seqnums)
        colnames.append('DAYOBS')
        units.append('None')
        columns.append(dayobs)
        colnames.append('FILE0')
        units.append('None')
        columns.append(ccd1.imfile.split('/')[-1])
        colnames.append('FILE1')
        units.append('None')
        columns.append(ccd2.imfile.split('/')[-1])
        formats = 'E'*(len(colnames) - 2) + 2*'J' + 2*'A'
        fits_cols = [fits.Column(name=colnames[i], format=formats[i],
                                 unit=units[i], array=columns[i])
                     for i in range(len(columns))]
        output.append(fitsTableFactory(fits_cols))
        output[-1].name = 'PTC_STATS'
        output[0].header['NAMPS'] = len(all_amps)
        fitsWriteto(output, outfile, overwrite=True)

    @staticmethod
    def fit_ptc_curve(mean, var, sig_cut=5, logger=None):
        """Fit the PTC curve for a set of mean-variance points."""
        index_old = []
        index = list(np.where((mean < 4e4)*(var >0))[0])
        count = 1
        # Initial guess for BF coeff, gain, and square of the read noise
        pars = 2.7e-6, 0.75, 25
        try:
            while index != index_old and count < 10:
                try:
                    results = scipy.optimize.leastsq(residuals, pars,
                                                     full_output=1,
                                                     args=(mean[index],
                                                           var[index]))
                except TypeError as err:
                    if logger is not None:
                        logger.info(err)
                        logger.info('Too few remaining mean-variance points:  %s' % len(index))

                pars, cov = results[:2]
                sig_resids = residuals(pars, mean, var)
                index_old = deepcopy(index)
                index = list(np.where(np.abs(sig_resids) < sig_cut)[0])
                count += 1

            ptc_a00 = pars[0]
            ptc_a00_error = np.sqrt(cov[0][0])
            ptc_gain = pars[1]
            ptc_error = np.sqrt(cov[1][1])
            ptc_noise = np.sqrt(pars[2])
            ptc_noise_error = 0.5/ptc_noise*np.sqrt(cov[2][2])
            # Cannot assume that the mean values are sorted
            ptc_turnoff = max(mean[index])
        except Exception as eobj:
            if logger is not None:
                logger.info("Exception caught while fitting PTC:")
                logger.info(str(eobj))
            ptc_gain = 0.
            ptc_error = -1.
            ptc_a00 = 0.
            ptc_a00_error = -1.
            ptc_noise = 0.
            ptc_noise_error = -1.
            ptc_turnoff = 0.
        return (ptc_gain, ptc_error, ptc_a00, ptc_a00_error, ptc_noise,
                ptc_noise_error, ptc_turnoff)

    def _fit_curves(self, ptc_stats, sensor_id, sig_cut=5):
        """
        Fit a model to the variance vs. mean
        """
        outfile = self.config.eotest_results_file
        if outfile is None:
            outfile = os.path.join(self.config.output_dir,
                                   '%s_eotest_results.fits' % sensor_id)
        print('Writing:  ' + outfile)
        output = EOTestResults(outfile, namps=len(ptc_stats))
        # Fit for gain and error, a00 and its uncertainty, inferred
        # noise and uncertainty, and the 'turnoff' level (in
        # electrons) and write to an EO test results file.
        for amp in ptc_stats:
            mean, var = np.array(ptc_stats[amp][0]), np.array(ptc_stats[amp][1])
            (ptc_gain, ptc_error, ptc_a00, ptc_a00_error, ptc_noise,
             ptc_noise_error, ptc_turnoff) \
             = self.fit_ptc_curve(mean, var, sig_cut=sig_cut, logger=self.log)
            output.add_seg_result(amp, 'PTC_GAIN', ptc_gain)
            output.add_seg_result(amp, 'PTC_GAIN_ERROR', ptc_error)
            output.add_seg_result(amp, 'PTC_A00', ptc_a00)
            output.add_seg_result(amp, 'PTC_A00_ERROR', ptc_a00_error)
            output.add_seg_result(amp, 'PTC_NOISE', ptc_noise)
            output.add_seg_result(amp, 'PTC_NOISE_ERROR', ptc_noise_error)
            output.add_seg_result(amp, 'PTC_TURNOFF', ptc_turnoff)
            self.log.info("%i  %f  %f %f %f %f %f %f" % (amp, ptc_gain,
                                                         ptc_error, ptc_a00,
                                                         ptc_a00_error,
                                                         ptc_noise,
                                                         ptc_noise_error,
                                                         ptc_turnoff))

            # Now fit the means and variances that were evaluated with MAD filtering
            mean, var = np.array(ptc_stats[amp][2]), np.array(ptc_stats[amp][3])
            (ptc_gain, ptc_error, ptc_a00, ptc_a00_error, ptc_noise,
             ptc_noise_error, ptc_turnoff) \
             = self.fit_ptc_curve(mean, var, sig_cut=sig_cut, logger=self.log)
            output.add_seg_result(amp, 'PTC_MGAIN', ptc_gain)
            output.add_seg_result(amp, 'PTC_MGAIN_ERROR', ptc_error)
            output.add_seg_result(amp, 'PTC_MA00', ptc_a00)
            output.add_seg_result(amp, 'PTC_MA00_ERROR', ptc_a00_error)
            output.add_seg_result(amp, 'PTC_MNOISE', ptc_noise)
            output.add_seg_result(amp, 'PTC_MNOISE_ERROR', ptc_noise_error)
            output.add_seg_result(amp, 'PTC_MTURNOFF', ptc_turnoff)
            self.log.info("%i  %f  %f %f %f %f %f %f" % (amp, ptc_gain,
                                                         ptc_error, ptc_a00,
                                                         ptc_a00_error,
                                                         ptc_noise,
                                                         ptc_noise_error,
                                                         ptc_turnoff))
        output.write()
