"""
brief Compute the brighter fatter correlation coefficients.  Input is
flat pairs.  Based on the LSST DM cdoe.  Output is the nearest neighbor
coefficients per amp.
"""
import os
import glob
from collections import namedtuple
import numpy as np
from astropy.io import fits
import lsst.geom as lsstGeom
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
    maxLag = pexConfig.Field("Maximum lag", int, default=2)
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
                          'COV10 COV10_err COV20 COV20_err COV01 COV01_err COV02 COV02_err COV11 COV11_err mean'.split())


class BFTask(pipeBase.Task):
    """Task to compute the brighter fatter correlation coefficients."""
    ConfigClass = BFConfig
    _DefaultName = "BrighterFatterTask"

    @pipeBase.timeMethod
    def run(self, sensor_id, flat_files, single_pairs=True,
            dark_frame=None, mask_files=(), flat2_finder=None,
            bias_frame=None, meanidx=0, linearity_correction=None,gains=None):
        """
        Compute the average nearest neighbor correlation coefficients for
        all flat pairs given for a particular exposure time.  Additionally
        store the flat means per amp.
        """

        if single_pairs:
            if flat2_finder is None:
                flat2_finder = find_flat2
            flats = find_flats(flat_files, flat2_finder=flat2_finder)
        else:
            flats = split_flats(flat_files)

        all_amps = imutils.allAmps(flat_files[0])

        # List with some number of flat pairs per exposure
        # [[(flat1,flat2),(flat1,flat2)],[(flat1,flat2),(flat1,flat2)]]

        BFResults = {amp: BFAmpResults([], [], [], [], [], [], [], [], [], [], []) for amp in all_amps}

        for flat_pair in flats:
            self.log.info("%s\n%s", *flat_pair)
            ccd1 = MaskedCCD(flat_pair[0], mask_files=mask_files,
                             bias_frame=bias_frame,
                             linearity_correction=linearity_correction,dark_frame=dark_frame)
            ccd2 = MaskedCCD(flat_pair[1], mask_files=mask_files,
                             bias_frame=bias_frame,
                             linearity_correction=linearity_correction,dark_frame=dark_frame)

            for amp in all_amps:
                self.log.info('on amp %s', amp)
                image1 = ccd1.unbiased_and_trimmed_image(amp)
                image2 = ccd2.unbiased_and_trimmed_image(amp)
                prepped_image1, mean1 = self.prep_image(image1, gains[amp])
                prepped_image2, mean2 = self.prep_image(image2, gains[amp])

                # Calculate the average mean of the pair.
                avemean = (mean1 + mean2)/2
                self.log.info('%s', avemean)

                # Compute the correlations.
                corr, corr_err = crossCorrelate(prepped_image1, prepped_image2,
                                                self.config.maxLag,
                                                self.config.nSigmaClip,
                                                self.config.backgroundBinSize)

                # Append the per-amp values for this pair to the
                # corresponding lists.
                BFResults[amp].COV10.append(corr[1][0])
                BFResults[amp].COV10_err.append(corr_err[1][0])
                BFResults[amp].COV20.append(corr[2][0])
                BFResults[amp].COV20_err.append(corr_err[2][0])
                BFResults[amp].COV01.append(corr[0][1])
                BFResults[amp].COV01_err.append(corr_err[0][1])
                BFResults[amp].COV02.append(corr[0][2])
                BFResults[amp].COV02_err.append(corr_err[0][2])
                BFResults[amp].COV11.append(corr[1][1])
                BFResults[amp].COV11_err.append(corr_err[1][1])
                BFResults[amp].mean.append(avemean)

        self.write_eotest_output(BFResults, sensor_id, meanidx=meanidx)

        return BFResults

    def fit_slopes(self,xcorr,ycorr,mean,adu_max):
        xcorr = np.array(xcorr)
        ycorr = np.array(ycorr)
        mean = np.array(mean)
        xcorr = xcorr[mean<adu_max]
        ycorr = ycorr[mean<adu_max]
        mean = mean[mean<adu_max]
        from scipy import stats
        slopex, _, _, _, errx = stats.linregress(mean, xcorr)
        slopey, _, _, _, erry = stats.linregress(mean, ycorr)
        return slopex, errx, slopey, erry

    def write_eotest_output(self, BFResults, sensor_id, meanidx=0, adu_max=1e5):
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
            colnames.extend(['AMP%02i_COV10' % amp, 'AMP%02i_COV10_ERR' % amp,
                             'AMP%02i_COV20' % amp, 'AMP%02i_COV20_ERR' % amp,
                             'AMP%02i_COV01' % amp, 'AMP%02i_COV01_ERR' % amp,
                             'AMP%02i_COV02' % amp, 'AMP%02i_COV02_ERR' % amp,
                             'AMP%02i_COV11' % amp, 'AMP%02i_COV11_ERR' % amp, 'AMP%02i_MEAN' % amp])
            units.extend(
                ['Unitless', 'Unitless', 'Unitless', 'Unitless', 'Unitless','Unitless','Unitless','Unitless','Unitless','Unitless','e-'])
            columns.extend([np.array(BFResults[amp][0], dtype=np.float),
                            np.array(BFResults[amp][1], dtype=np.float),
                            np.array(BFResults[amp][2], dtype=np.float),
                            np.array(BFResults[amp][3], dtype=np.float),
                            np.array(BFResults[amp][4], dtype=np.float),
                            np.array(BFResults[amp][5], dtype=np.float),
                            np.array(BFResults[amp][6], dtype=np.float),
                            np.array(BFResults[amp][7], dtype=np.float),
                            np.array(BFResults[amp][8], dtype=np.float),
                            np.array(BFResults[amp][9], dtype=np.float),
                            np.array(BFResults[amp][10], dtype=np.float)])
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
            slopex, slopex_err, slopey_err, slopey \
                = self.fit_slopes(BFResults[amp][0], BFResults[amp][2],
                                  BFResults[amp][4], adu_max)
            results.add_seg_result(amp, 'BF_SLOPEX', slopex)
            results.add_seg_result(amp, 'BF_SLOPEX_ERR', slopex_err)
            results.add_seg_result(amp, 'BF_SLOPEY', slopey)
            results.add_seg_result(amp, 'BF_SLOPEY_ERR', slopey_err)

        results.write(clobber=True)

    def prep_image(self, exp, gain):
        """
        Crop the image to avoid edge effects based on the Config border
        parameter. Additionally, if there is a dark image, subtract.
        """

        # Clone so that we don't modify the input image.
        local_exp = exp.clone()

        # The border that we wish to crop.
        border = self.config.nPixBorder

        sctrl = afwMath.StatisticsControl()

        # Crop the image within a border region.
        bbox = local_exp.getBBox()
        bbox.grow(-border)
        local_exp = local_exp[bbox]
        local_exp*=gain

        # Calculate the mean of the image.
        mean = afwMath.makeStatistics(local_exp, afwMath.MEDIAN,
                                      sctrl).getValue()

        return local_exp, mean


def crossCorrelate(maskedimage1, maskedimage2, maxLag, sigma, binsize):
    """
    Calculate the correlation coefficients.
    """
    sctrl = afwMath.StatisticsControl()
    sctrl.setNumSigmaClip(sigma)
    mask = maskedimage1.getMask()
    INTRP = mask.getPlaneBitMask("INTRP")
    sctrl.setAndMask(INTRP)


    # Diff the images.
    diff = maskedimage1.clone()
    diff -= maskedimage2.getImage()

    # Subtract background.
    nx = diff.getWidth()//binsize
    ny = diff.getHeight()//binsize
    bctrl = afwMath.BackgroundControl(nx, ny, sctrl, afwMath.MEDIAN)
    bkgd = afwMath.makeBackground(diff, bctrl)
    bgImg = bkgd.getImageF(afwMath.Interpolate.CUBIC_SPLINE,
                           afwMath.REDUCE_INTERP_ORDER)

    diff -= bgImg

    # Measure the correlations
    x0, y0 = diff.getXY0()
    width, height = diff.getDimensions()
    bbox_extent = lsstGeom.Extent2I(width - maxLag, height - maxLag)

    bbox = lsstGeom.Box2I(lsstGeom.Point2I(x0, y0), bbox_extent)
    dim0 = diff[bbox].clone()
    dim0 -= afwMath.makeStatistics(dim0, afwMath.MEDIAN, sctrl).getValue()

    xcorr = np.zeros((maxLag + 1, maxLag + 1), dtype=np.float64)
    xcorr_err = np.zeros((maxLag + 1, maxLag + 1), dtype=np.float64)

    for xlag in range(maxLag + 1):
        for ylag in range(maxLag + 1):
            bbox_lag = lsstGeom.Box2I(lsstGeom.Point2I(x0 + xlag, y0 + ylag),
                                      bbox_extent)
            dim_xy = diff[bbox_lag].clone()
            dim_xy -= afwMath.makeStatistics(dim_xy, afwMath.MEDIAN,
                                             sctrl).getValue()
            dim_xy *= dim0
            xcorr[xlag, ylag] = afwMath.makeStatistics(
                dim_xy, afwMath.MEDIAN, sctrl).getValue()
            dim_xy_array = dim_xy.getImage().getArray().flatten()/xcorr[0][0]
            N = len(dim_xy_array.flatten())
            if xlag != 0 and ylag != 0:
                f = (1+xcorr[xlag, ylag]/xcorr[0][0]) / \
                    (1-xcorr[xlag, ylag]/xcorr[0][0])
                xcorr_err[xlag, ylag] = (
                    np.std(dim_xy_array)/np.sqrt(N))*np.sqrt(f)
            else:
                xcorr_err[xlag, ylag] = 0
    return xcorr, xcorr_err
