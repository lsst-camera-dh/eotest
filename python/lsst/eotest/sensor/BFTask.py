"""
@brief Compute the brighter fatter correlation coefficients.  Input is
flat pairs.  Output is the nearest neighbor coefficients per amp.

@author E. Phillips Longley <elp25@duke.edu>
"""
import os
import glob
import numpy as np

import astropy.io.fits as fits
from lsst.eotest.fitsTools import fitsTableFactory, fitsWriteto
import lsst.eotest.image_utils as imutils
from lsst.eotest.sensor.MaskedCCD import MaskedCCD

import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

import lsst.eotest.image_utils as im_util
from lsst.eotest.sensor import MaskedCCD, AmplifierGeometry
from lsst.eotest.sensor.AmplifierGeometry import makeAmplifierGeometry
from lsst.eotest.sensor.EOTestResults import EOTestResults

def find_flat2(flat1):
    """Given flat1, finds flat2."""
    pattern = flat1.split('flat1')[0] + 'flat2*.fits'
    try:
        flat2 = glob.glob(pattern)[0]
        return flat2
    except IndexError:
        return flat1

def split_flats(flats):
    """Given a set of flats, splits them into flat1s and flat2s consecutively.
    """
    flat1s = [flats[i] for i in range(len(flats)) if i%2==0]
    flat2s = [flats[i] for i in range(len(flats)) if i%2==1]
    return [(flat1s[i],flat2s[i]) for i in range(len(flat1s))]

def glob_flats(full_path, outfile='ptc_flats.txt'):
    flats = glob.glob(os.path.join(full_path, '*_flat?.fits'))
    output = open(outfile, 'w')
    for item in flats:
        output.write('%s\n' % item)
    output.close()

def find_flats(flats):
    file1s = sorted([item.strip() for item in flats
                     if item.find('flat1') != -1])
    return [[(f1,find_flat2(f1))] for f1 in file1s]

class BFConfig(pexConfig.Config):
    """Configuration for BFTask"""
    maxLag = pexConfig.Field("Maximum lag",
                             int, default=8)
    nPixBorder = pexConfig.Field("Number of pixels to clip on the border",
                                int, default=10)
    nSigmaClip = pexConfig.Field("Number of sigma to clip for corr calc", int,
    default = 3)
    output_dir = pexConfig.Field("Output directory", str, default="/u/ec/elp25/private/TESTBFOUTPUT/")
    eotest_results_file = pexConfig.Field("EO test results filename",
                                          str, default=None)
    verbose = pexConfig.Field("Turn verbosity on", bool, default=True)
    backgroundBinSize = pexConfig.Field("Background bin size",int, default = 128)


class BFResults(object):
    def __init__(self, C10s, C01s, means):
        self.C10s = C10s
        self.C01s = C01s
        self.means = means

class BFTask(pipeBase.Task):
    """Task to compute the brighter fatter correlation coefficients."""
    ConfigClass = BFConfig
    _DefaultName = "BrighterFatterTask"

    @pipeBase.timeMethod

    def run(self, sensor_id, flat_files,meanidx=0,single_pairs=True, bias_frames=None,dark_frame = None, mask_files = ()):
        """Compute the average nearest neighbor correlation coefficients for
        all flat pairs given for a particular exposure time.  Additionally
        store the flat means per amp."""
        if single_pairs:
            flats = find_flats(flat_files)
        else:
            flats = split_flats(flat_files)

        all_amps = imutils.allAmps(flat_files[0])

        # List with some number of flat pairs per exposure
        # [[(flat1,flat2),(flat1,flat2)],[(flat1,flat2),(flat1,flat2)]]

        C10s = {}
        C01s = {}
        means = {}
        BFResults = {}

        for amp in all_amps:
            print('on amp',amp)
            C10_exp = []
            C01_exp = []
            mean_exp = []

            for exposure in range(len(flats)):
                C10_amp = []
                C01_amp = []
                mean_amp = []
                if dark_frame:
                    ccd_dark = MaskedCCD(dark_frame,mask_files = mask_files)
                    dark_frame = ccd_dark.unbiased_and_trimmed_image(amp)
                for flat_pair in flats[exposure]:
                    print(flat_pair)
                    ccd1 = MaskedCCD(flat_pair[0],mask_files = mask_files)
                    ccd2 = MaskedCCD(flat_pair[1],mask_files = mask_files)

                    image1 = ccd1.unbiased_and_trimmed_image(amp)
                    image2 = ccd2.unbiased_and_trimmed_image(amp)
                    prepped_image1, mean1 = self.prep_image(image1,dark_frame)
                    prepped_image2, mean2 = self.prep_image(image2,dark_frame)
                    #Calculate the average mean of the pair
                    avemean = (mean1+mean2)/2

                    #Compute the correlations
                    corr = self.crossCorrelate(prepped_image1,prepped_image2)
                    C10_amp.append(corr[1][0]/corr[0][0])
                    C01_amp.append(corr[0][1]/corr[0][0])
                    mean_amp.append(avemean)

                C10_exp.append(np.mean(C10_amp))
                C01_exp.append(np.mean(C01_amp))
                mean_exp.append(np.mean(mean_amp))

            #Store cov/var
            C10s[amp] = C10_exp
            C01s[amp] = C01_exp
            means[amp] = mean_exp
            BFResults[amp] = [C10_exp,C01_exp,mean_exp]

        self.writeEotestOutput(BFResults, all_amps,sensor_id,meanidx)

        return


    def writeEotestOutput(self,BFResults,all_amps,sensor_id,meanidx=0):

        outfile = os.path.join(self.config.output_dir,
                               '%s_bftest.fits' % sensor_id)
        output = fits.HDUList()
        output.append(fits.PrimaryHDU())
        colnames = []
        units = []
        columns = []
        for amp in all_amps:
            colnames.extend(['AMP%02i_C10' % amp, 'AMP%02i_C01' % amp,
            'AMP%02i_MEAN' % amp])
            units.extend(['Unitless', 'Unitless','Unitless'])
            columns.extend([np.array(BFResults[amp][0], dtype=np.float),
                            np.array(BFResults[amp][1], dtype=np.float),
                            np.array(BFResults[amp][2], dtype=np.float)])
        formats = 'E'*len(colnames)
        fits_cols = [fits.Column(name=colnames[i], format=formats[i],
                                 unit=units[i], array=columns[i])
                     for i in range(len(columns))]
        output.append(fitsTableFactory(fits_cols))
        output[-1].name = 'BF_STATS'
        output[0].header['NAMPS'] = len(all_amps)
        fitsWriteto(output, outfile, clobber=True)
        #Output a file of the coefficients at a given mean, given
        #as the index of the exposure in the list.

        results_file = os.path.join(self.config.output_dir,
                                        '%s_eotest_results.fits' % sensor_id)


        results = EOTestResults(results_file, namps=len(all_amps))

        for amp in all_amps:
            results.add_seg_result(amp, 'C10', BFResults[amp][0][meanidx])
            results.add_seg_result(amp, 'C01', BFResults[amp][1][meanidx])
            results.add_seg_result(amp, 'MEAN', BFResults[amp][2][meanidx])


        results.write(clobber=True)

    def prep_image(self,exp,dark_image):
        """Crop the image to avoid edge effects based on the Config border
        parameter. Additionally, if there is a dark image, subtract."""

        # Clone so that we don't modify the input image.
        local_exp = exp.clone()
        del exp

        # The border that we wish to crop.
        border = self.config.nPixBorder

        sctrl = afwMath.StatisticsControl()

        # If a dark image is passed, subtract it.
        if dark_image:
            local_exp -= dark_image.getImage()

        # Calculate the mean of the image.
        mean = afwMath.makeStatistics(local_exp[border:-border, border:-border],
                                                            afwMath.MEAN, sctrl).getValue()
        # Crop the image.
        local_exp = local_exp[border:-border, border:-border]

        return local_exp, mean

    def crossCorrelate(self,maskedimage1, maskedimage2):

        maxLag = self.config.maxLag
        border = self.config.nPixBorder
        sigma = self.config.nSigmaClip

        sctrl = afwMath.StatisticsControl()
        sctrl.setNumSigmaClip(sigma)

        # Diff the images, and apply border
        diff = maskedimage1.clone()
        diff -= maskedimage2.getImage()

        # Subtract background.  It should be a constant, but it isn't always
        binsize = self.config.backgroundBinSize
        nx = diff.getWidth()//binsize
        ny = diff.getHeight()//binsize
        bctrl = afwMath.BackgroundControl(nx, ny, sctrl, afwMath.MEANCLIP)
        bkgd = afwMath.makeBackground(diff, bctrl)
        bgImg = bkgd.getImageF(afwMath.Interpolate.CUBIC_SPLINE, afwMath.REDUCE_INTERP_ORDER)
        bgMean = np.mean(bgImg.getArray())

        diff -= bgImg

        # Measure the correlations
        dim0 = diff[0: -maxLag, : -maxLag]
        dim0 -= afwMath.makeStatistics(dim0, afwMath.MEANCLIP, sctrl).getValue()
        width, height = dim0.getDimensions()
        xcorr = np.zeros((maxLag + 1, maxLag + 1), dtype=np.float64)

        for xlag in range(maxLag + 1):
            for ylag in range(maxLag + 1):
                dim_xy = diff[xlag:xlag + width, ylag: ylag + height].clone()
                dim_xy -= afwMath.makeStatistics(dim_xy, afwMath.MEANCLIP, sctrl).getValue()
                dim_xy *= dim0
                xcorr[xlag, ylag] = afwMath.makeStatistics(dim_xy,afwMath.MEANCLIP, sctrl).getValue()

        return xcorr
