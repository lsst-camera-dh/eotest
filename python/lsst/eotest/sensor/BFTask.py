"""
brief Compute the brighter fatter correlation coefficients.  Input is
flat pairs.  Based on the LSST DM cdoe.  Output is the nearest neighbor
coefficients per amp.
"""
import os
import glob
from collections import namedtuple
import subprocess
import psutil
import numpy as np
from astropy.io import fits
import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.eotest import fitsTools
import lsst.eotest.image_utils as imutils
from .MaskedCCD import MaskedCCD
from .EOTestResults import EOTestResults


def find_flat2(flat1):
    """Given flat1, finds flat2."""
    pattern = flat1.split('flat1')[0] + 'flat2*.fits'
    try:
        flat2 = glob.glob(pattern)[0]
        return flat2
    except IndexError:
        return flat1


def split_flats(flats):
    """
    Given a set of flats, splits them into flat1s and flat2s consecutively.
    """
    flat1s = [flats[i] for i in range(len(flats)) if i % 2 == 0]
    flat2s = [flats[i] for i in range(len(flats)) if i % 2 == 1]
    return [(flat1s[i], flat2s[i]) for i in range(len(flat1s))]


def find_flats(flats, flat2_finder=find_flat2):
    """Find flat pairs."""
    file1s = sorted([item.strip() for item in flats
                     if item.find('flat1') != -1])
    return [(f1, flat2_finder(f1)) for f1 in file1s]


class BFConfig(pexConfig.Config):
    """Configuration for BFTask"""
    maxLag = pexConfig.Field("Maximum lag", int, default=1)
    nPixBorder = pexConfig.Field("Number of pixels to clip on the border",
                                 int, default=10)
    nSigmaClip = pexConfig.Field("Number of sigma to clip for corr calc", int,
                                 default=3)
    output_dir = pexConfig.Field("Output directory", str, default=".")
    eotest_results_file = pexConfig.Field("EO test results filename", str,
                                          default=None)
    verbose = pexConfig.Field("Turn verbosity on", bool, default=True)
    backgroundBinSize = pexConfig.Field("Background bin size", int, default=128)


BFAmpResults = namedtuple('BFAmpResults',
                          'xcorr xcorr_err ycorr ycorr_err mean'.split())


class BFTask(pipeBase.Task):
    """Task to compute the brighter fatter correlation coefficients."""
    ConfigClass = BFConfig
    _DefaultName = "BrighterFatterTask"

    def log_mem_info(self, message, amp=None):
        if self.process is None:
            return
        rss_mem = self.process.memory_info().rss/1024.**3
        amp_info = ", amp %d" % amp if amp is not None else ''
        self.log.info('BFTask, %s: rss_mem=%.2f GB; %s',
                      self.process.pid, rss_mem, message + amp_info)

    def _set_process(self):
        if os.environ.get('LCATR_LOG_MEM_INFO', False) == 'True':
            self.process = psutil.Process(os.getpid())
        else:
            self.process = None

    @pipeBase.timeMethod
    def run(self, sensor_id, flat_files, single_pairs=True,
            dark_frame=None, mask_files=(), flat2_finder=None,
            bias_frame=None, meanidx=0):
        """
        Compute the average nearest neighbor correlation coefficients for
        all flat pairs given for a particular exposure time.  Additionally
        store the flat means per amp.
        """
        self._set_process()
        self.log_mem_info('entered run method')
        command = '(ulimit -a) >& ulimit_info_{}.txt'.format(sensor_id)
        subprocess.check_call(command, shell=True)
        if dark_frame is not None:
            self.log_mem_info('creating ccd_dark')
            ccd_dark = MaskedCCD(dark_frame, mask_files=mask_files)

        if single_pairs:
            if flat2_finder is None:
                flat2_finder = find_flat2
            self.log_mem_info('calling find_flats')
            flats = find_flats(flat_files, flat2_finder=flat2_finder)
        else:
            self.log_mem_info('calling split_flats')
            flats = split_flats(flat_files)

        self.log_mem_info('calling imutils.allAmps')
        all_amps = imutils.allAmps(flat_files[0])

        # List with some number of flat pairs per exposure
        # [[(flat1,flat2),(flat1,flat2)],[(flat1,flat2),(flat1,flat2)]]

        BFResults = {amp: BFAmpResults([], [], [], [], []) for amp in all_amps}

        for flat_pair in flats:
            self.log.info("%s\n%s", *flat_pair)
            self.log_mem_info('making ccd1')
            ccd1 = MaskedCCD(flat_pair[0], mask_files=mask_files,
                             bias_frame=bias_frame)
            self.log_mem_info('making ccd2')
            ccd2 = MaskedCCD(flat_pair[1], mask_files=mask_files,
                             bias_frame=bias_frame)

            for amp in all_amps:
                self.log.info('on amp %s', amp)
                self.log_mem_info('making dark_image', amp=amp)
                dark_image = None if dark_frame is None \
                             else ccd_dark.unbiased_and_trimmed_image(amp)

                self.log_mem_info('making image1', amp=amp)
                image1 = ccd1.unbiased_and_trimmed_image(amp)
                self.log_mem_info('making image2', amp=amp)
                image2 = ccd2.unbiased_and_trimmed_image(amp)
                self.log_mem_info('calling prep_image(image1)', amp=amp)
                prepped_image1, mean1 = self.prep_image(image1, dark_image)
                self.log_mem_info('calling prep_image(image2)', amp=amp)
                prepped_image2, mean2 = self.prep_image(image2, dark_image)

                # Calculate the average mean of the pair.
                avemean = (mean1 + mean2)/2
                self.log.info('%s', avemean)

                # Compute the correlations.
                self.log_mem_info('calling crossCorrelate', amp=amp)
                corr, corr_err \
                    = self.crossCorrelate(prepped_image1, prepped_image2,
                                          self.config.maxLag,
                                          self.config.nSigmaClip,
                                          self.config.backgroundBinSize)

                # Append the per-amp values for this pair to the
                # corresponding lists.
                BFResults[amp].xcorr.append(corr[1][0]/corr[0][0])
                BFResults[amp].xcorr_err.append(corr_err[1][0])
                BFResults[amp].ycorr.append(corr[0][1]/corr[0][0])
                BFResults[amp].ycorr_err.append(corr_err[0][1])
                BFResults[amp].mean.append(avemean)

        self.log_mem_info('calling write_eotest_output')
        self.write_eotest_output(BFResults, sensor_id, meanidx=meanidx)

        return BFResults

    def write_eotest_output(self, BFResults, sensor_id, meanidx=0):
        """Write the correlation curves to a FITS file for plotting,
        and the BF results to the eotest results file."""
        outfile = os.path.join(self.config.output_dir,
                               '%s_bf.fits' % sensor_id)
        output = fits.HDUList()
        output.append(fits.PrimaryHDU())
        colnames = []
        units = []
        columns = []
        for amp in BFResults:
            colnames.extend(['AMP%02i_XCORR' % amp, 'AMP%02i_XCORR_ERR' % amp,
                             'AMP%02i_YCORR' % amp, 'AMP%02i_YCORR_ERR' % amp,
                             'AMP%02i_MEAN' % amp])
            units.extend(
                ['Unitless', 'Unitless', 'Unitless', 'Unitless', 'ADU'])
            columns.extend([np.array(BFResults[amp][0], dtype=np.float),
                            np.array(BFResults[amp][1], dtype=np.float),
                            np.array(BFResults[amp][2], dtype=np.float),
                            np.array(BFResults[amp][3], dtype=np.float),
                            np.array(BFResults[amp][4], dtype=np.float)])
        formats = 'E'*len(colnames)
        fits_cols = [fits.Column(name=colnames[i], format=formats[i],
                                 unit=units[i], array=columns[i])
                     for i in range(len(columns))]
        output.append(fitsTools.fitsTableFactory(fits_cols))
        output[-1].name = 'BF_STATS'
        output[0].header['NAMPS'] = len(BFResults)
        fitsTools.fitsWriteto(output, outfile, clobber=True)

        # Output a file of the coefficients at a given mean, given
        # as the index of the exposure in the list.
        results_file = self.config.eotest_results_file
        if results_file is None:
            results_file = os.path.join(self.config.output_dir,
                                        '%s_eotest_results.fits' % sensor_id)

        results = EOTestResults(results_file, namps=len(BFResults))

        for amp in BFResults:
            results.add_seg_result(amp, 'BF_XCORR', BFResults[amp][0][meanidx])
            results.add_seg_result(amp, 'BF_XCORR_ERR', BFResults[amp][1][meanidx])
            results.add_seg_result(amp, 'BF_YCORR', BFResults[amp][2][meanidx])
            results.add_seg_result(amp, 'BF_YCORR_ERR', BFResults[amp][3][meanidx])
            results.add_seg_result(amp, 'BF_MEAN', BFResults[amp][4][meanidx])

        results.write(clobber=True)

    def prep_image(self, exp, dark_image=None):
        """
        Crop the image to avoid edge effects based on the Config border
        parameter. Additionally, if there is a dark image, subtract.
        """

        # Clone so that we don't modify the input image.
        local_exp = exp.clone()

        # The border that we wish to crop.
        border = self.config.nPixBorder

        sctrl = afwMath.StatisticsControl()

        # If a dark image is passed, subtract it.
        if dark_image is not None:
            local_exp -= dark_image.getImage()

        # Crop the image within a border region.
        bbox = local_exp.getBBox()
        bbox.grow(-border)
        local_exp = local_exp[bbox]

        # Calculate the mean of the image.
        mean = afwMath.makeStatistics(local_exp, afwMath.MEAN,
                                      sctrl).getValue()

        return local_exp, mean

    def crossCorrelate(self, maskedimage1, maskedimage2, maxLag, sigma,
                       binsize):
        """
        Calculate the correlation coefficients.
        """
        process = psutil.Process(os.getpid())
        self.log_mem_info('entered crossCorrelate')
        sctrl = afwMath.StatisticsControl()
        sctrl.setNumSigmaClip(sigma)

        # Diff the images.
        self.log_mem_info('maskedimage1.clone()')
        diff = maskedimage1.clone()
        self.log_mem_info('diff -= maskedimage2.getImage()')
        diff -= maskedimage2.getImage()

        # Subtract background.
        nx = diff.getWidth()//binsize
        ny = diff.getHeight()//binsize
        self.log_mem_info('creating bctrl')
        bctrl = afwMath.BackgroundControl(nx, ny, sctrl, afwMath.MEANCLIP)
        self.log_mem_info('creating bkgd')
        bkgd = afwMath.makeBackground(diff, bctrl)
        self.log_mem_info('creating bgImg')
        bgImg = bkgd.getImageF(afwMath.Interpolate.CUBIC_SPLINE,
                               afwMath.REDUCE_INTERP_ORDER)

        self.log_mem_info('diff -= bgImg')
        diff -= bgImg

        # Measure the correlations
        x0, y0 = diff.getXY0()
        width, height = diff.getDimensions()
        bbox_extent = afwGeom.Extent2I(width - maxLag, height - maxLag)

        bbox = afwGeom.Box2I(afwGeom.Point2I(x0, y0), bbox_extent)
        dim0 = diff[bbox].clone()
        self.log_mem_info('dim0 -= afwMath.makeStatistics(...)')
        dim0 -= afwMath.makeStatistics(dim0, afwMath.MEANCLIP,
                                       sctrl).getValue()

        xcorr = np.zeros((maxLag + 1, maxLag + 1), dtype=np.float64)
        xcorr_err = np.zeros((maxLag + 1, maxLag + 1), dtype=np.float64)

        for xlag in range(maxLag + 1):
            for ylag in range(maxLag + 1):
                self.log_mem_info('xlag = %d, ylag = %d' % (xlag, ylag))
                bbox_lag = afwGeom.Box2I(afwGeom.Point2I(x0 + xlag, y0 + ylag),
                                         bbox_extent)
                dim_xy = diff[bbox_lag].clone()
                self.log_mem_info('dim_xy -= afwMath.makeStatistics(...)')
                dim_xy -= afwMath.makeStatistics(dim_xy, afwMath.MEANCLIP,
                                                 sctrl).getValue()
                dim_xy *= dim0
                self.log_mem_info('computing xcorr[xlag, ylag]')
                xcorr[xlag, ylag] = afwMath.makeStatistics(
                    dim_xy, afwMath.MEANCLIP, sctrl).getValue()
                self.log_mem_info('dim_xy_image =')
                dim_xy_image = dim_xy.getImage()
                self.log_mem_info('dim_xy_array =')
                dim_xy_array = dim_xy_image.array/xcorr[0][0]
                self.log_mem_info('N = np.prod(...)')
                N = np.prod(dim_xy_array.shape)
                if xlag != 0 and ylag != 0:
                    self.log_mem_info('f = (1+xcorr...)')
                    f = (1+xcorr[xlag, ylag]/xcorr[0][0]) / \
                        (1-xcorr[xlag, ylag]/xcorr[0][0])
                    xcorr_err[xlag, ylag] \
                        = np.std(dim_xy_array)/np.sqrt(N)*np.sqrt(f)
                else:
                    xcorr_err[xlag, ylag] = 0
        return xcorr, xcorr_err
