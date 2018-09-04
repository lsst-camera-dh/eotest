import numpy as np
import warnings
import itertools
from astropy.io import fits
from lsst.eotest.fitsTools import fitsWriteto, fitsTableFactory

import lsst.afw.detection as afwDetect
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath

import scipy.optimize
from scipy.special import erf, gammaincc

import lsst.eotest.image_utils as imutils
from MaskedCCD import MaskedCCD, MaskedCCDBiasImageException
from fe55_psf import prect_values, chisq, psf_func, residuals, pixel_integral
from AmplifierGeometry import parse_geom_kwd

_sqrt2 = np.sqrt(2)

def get_amp_number(amp_coords, x, y):
    """Determine amplifier number from x, y mosaic coordinates."""

    for key in amp_coords.keys():
        x_lims = amp_coords[key][:2]
        y_lims = amp_coords[key][2:]

        if min(x_lims) <= x < max(x_lims) and min(y_lims) <= y < max(y_lims):
            return key

class SpotMomentFit2(object):
    """Copy of spot moment fitter for use on combined CCD image."""

    def __init__(self, nsig=10, outfile=None):

        self.nsig = nsig
        self.outfile = outfile
        if outfile is None:
            self.output = fits.HDUList()
            self.output.append(fits.PrimaryHDU())
        else:
            self.output = fits.open(self.outfile)

    def process_image(self, ccd, infile, sigma0=0.36, dn0=1590./5.,
                      bg_reg=(10, 10), logger=None, oscan_fit_order=1):

        foo = fits.open(infile)
        datasec = parse_geom_kwd(foo[1].header['DATASEC'])
        nx_segments = 8
        ny_segments = len(ccd)/nx_segments
        nx = nx_segments*(datasec['xmax'] - datasec['xmin'] + 1)
        ny = ny_segments*(datasec['ymax'] - datasec['ymin'] + 1)
        mosaic = np.zeros((ny, nx), dtype=np.float32)
        amp_coords = {}
        for ypos in range(ny_segments):
            for xpos in range(nx_segments):
                amp = ypos*nx_segments + xpos + 1
                #
                # Determine subarray boundaries in the mosaicked image array
                # from DETSEC keywords for each segment.
                detsec = parse_geom_kwd(foo[amp].header['DETSEC'])
                xmin = nx - max(detsec['xmin'], detsec['xmax'])
                xmax = nx - min(detsec['xmin'], detsec['xmax']) + 1
                ymin = ny - max(detsec['ymin'], detsec['ymax'])
                ymax = ny - min(detsec['ymin'], detsec['ymax']) + 1
                #
                # Save coordinates of segment for later use
                #
                amp_coords[amp] = xmin, xmax, ymin, ymax
                #
                # Extract the bias-subtracted masked image for this segment.
                segment_image = ccd.unbiased_and_trimmed_image(amp)
                subarr = segment_image.getImage().getArray()
                #
                # Determine flips in x- and y-directions in order to
                # get the (1, 1) pixel in the lower right corner.
                if detsec['xmax'] > detsec['xmin']:  # flip in x-direction
                    subarr = subarr[:, ::-1]
                if detsec['ymax'] > detsec['ymin']:  # flip in y-direction
                    subarr = subarr[::-1, :]
                #
                # Set the subarray in the mosaicked image.
                mosaic[ymin:ymax, xmin:xmax] = subarr

        for key in amp_coords.keys():
            print key, amp_coords[key]

        image = afwImage.ImageF(mosaic)
        imarr = image.getArray()

        flags = afwMath.MEDIAN | afwMath.STDEVCLIP
        statistics = afwMath.makeStatistics(image, flags, ccd.stat_ctrl)
        median = statistics.getValue(afwMath.MEDIAN)
        stdev = statistics.getValue(afwMath.STDEVCLIP)

        threshold = afwDetect.Threshold(median + self.nsig*stdev)
        npix_min = 60
        if logger is not None:
            logger.info("SpotMomentFit.process_image: threshold= %s"
                        % threshold.getValue())
        fp_set = afwDetect.FootprintSet(image, threshold, npix_min)

        ## Iterate over each footprint and calculate moments
        amp = [] # Array for amplifier location identifier
        x0, y0 = [], [] # These are arrays to hold per-spot measurements
        prect_data = [] # This holds flattened spot footprint
        sigmax, sigmay, dn, dn_fp, chiprob = [], [], [], [], [] # Arrays to hold the moments fitting measurements
        chi2s, dofs = [], []
        maxDNs = []
        failed_curve_fits = 0
        num_fp = 0
        for fp in fp_set.getFootprints():

            num_fp +=1
            peak = [pk for pk in fp.getPeaks()][0]
            spans = fp.getSpans()
            positions = []
            zvals = []
            dn_sum = 0

            for span in spans:
                y = span.getY()
                for x in range(span.getX0(), span.getX1() + 1):
                    positions.append((x, y))
                    zvals.append(imarr[y][x])
                    dn_sum += imarr[y][x]

            try:
                # Use clipped stdev as DN error estimate for all pixels
                dn_errors = stdev*np.ones(len(positions))
                p0 = (peak.getIx(), peak.getIy(), sigma0, sigma0, dn0)
                pars, _ = scipy.optimize.leastsq(residuals, p0,
                                                     args=(positions, zvals,
                                                           dn_errors))
                sigmax.append(np.abs(pars[2]))
                sigmay.append(np.abs(pars[3]))
                dn.append(pars[4])
                x0.append(pars[0])
                y0.append(pars[1])
                amp.append(get_amp_number(amp_coords, pars[0], pars[1]))

                ## Determine amplifier position
                fp_amp = get_amp_number(amp_coords, x0, y0)

                dn_fp.append(dn_sum)
                chi2 = chisq(positions, zvals, x0[-1], y0[-1],
                             sigmax[-1], sigmay[-1], dn[-1], dn_errors)
                dof = fp.getNpix() - 5 # number of pars for our fit
                chiprob.append(gammaincc(dof/2., chi2/2.))
                chi2s.append(chi2)
                dofs.append(dof)
                maxDNs.append(max(zvals))

                prect_data_row = prect_values(peak,imarr, ixm=10, ixp=10, iym=10, iyp=10)
                prect_data.append(prect_data_row)
            except RuntimeError:
                failed_curve_fits +=1

        self._save_ext_data(amp, x0, y0,prect_data, sigmax, sigmay, dn, dn_fp, chiprob,
                            chi2s, dofs, maxDNs)

        print(num_fp)

    def _save_ext_data(self, amp, x0, y0, prect_data, sigmax, sigmay, dn, dn_fp, chiprob,
                        chi2s, dofs, maxDNs):
        """
        Write results from the source detection and Gaussian fitting
        to the FITS extension corresponding to the specified
        amplifier.
        """
        extname = 'MOSAIC'
        try:
            #
            # Append new rows if HDU for this segment already exists.
            #
            table_hdu = self.output[extname]
            row0 = table_hdu.header['NAXIS2']
            nrows = row0 + len(x0)
            table_hdu = fitsTableFactory(table_hdu.data, nrows=nrows)
            for i in range(len(x0)):
                row = i + row0
                table_hdu.data[row]['AMPLIFIER'] = amp[i]
                table_hdu.data[row]['XPOS'] = x0[i]
                table_hdu.data[row]['YPOS'] = y0[i]
                table_hdu.data[row]['PRECT_DATA'] = prect_data[i]
                table_hdu.data[row]['SIGMAX'] = sigmax[i]
                table_hdu.data[row]['SIGMAY'] = sigmay[i]
                table_hdu.data[row]['DN'] = dn[i]
                table_hdu.data[row]['DN_FP_SUM'] = dn_fp[i]
                table_hdu.data[row]['CHIPROB'] = chiprob[i]
                table_hdu.data[row]['CHI2'] = chi2s[i]
                table_hdu.data[row]['DOF'] = dofs[i]
                table_hdu.data[row]['MAXDN'] = maxDNs[i]
            table_hdu.name = extname
            self.output[extname] = table_hdu
        except KeyError:
            #
            # Extension for this segment does not yet exist, so add it.
            #
            colnames = ['AMPLIFIER', 'XPOS', 'YPOS', 'PRECT_DATA','SIGMAX', 'SIGMAY', 'DN',
                        'DN_FP_SUM', 'CHIPROB', 'CHI2', 'DOF', 'MAXDN']

            columns = [np.array(amp), np.array(x0), np.array(y0), np.array(prect_data),
                       np.array(sigmax), np.array(sigmay),
                       np.array(dn), np.array(dn_fp),
                       np.array(chiprob), np.array(chi2s), np.array(dofs),
                       np.array(maxDNs)]
            formats  = ['I'] + ['E']*(2) + ['441E'] +['E']*(8)
            units = ['None', 'pixel', 'pixel', 'pixel', 'pixel',
                     'ADU', 'ADU', 'None', 'None', 'None', 'ADU',
                     'pixel', 'pixel', 'ADU']
            fits_cols = lambda coldata: [fits.Column(name=colname,
                                                     format=format,
                                                     unit=unit,
                                                     array=column)
                                         for colname, format, unit, column
                                         in coldata]
            self.output.append(fitsTableFactory(fits_cols(zip(colnames,
                                                              formats,
                                                              units,
                                                              columns))))
            self.output[-1].name = extname

    def write_results(self, outfile='spot_moment_params.fits'):
        fitsWriteto(self.output, outfile, overwrite=True)

class SpotMomentFit(object):

    def __init__(self, nsig=10, outfile=None):

        self.nsig = nsig

        self.amp_set = set()
        self.outfile = outfile
        if outfile is None:
            self.output = fits.HDUList()
            self.output.append(fits.PrimaryHDU())
        else:
            # Append new data to existing file.
            self.output = fits.open(self.outfile)

    def _bg_image(self, ccd, amp, nx, ny):
        "Compute background image based on clipped local mean."
        bg_ctrl = afwMath.BackgroundControl(nx, ny, ccd.stat_ctrl)
        bg = afwMath.makeBackground(ccd[amp], bg_ctrl)
        return bg.getImageF()

    def process_image(self, ccd, amp, sigma0=0.36, dn0=1590./5.,
                      bg_reg=(10, 10), logger=None, oscan_fit_order=1):

        try:
            image = ccd.bias_subtracted_image(amp, fit_order=oscan_fit_order)
        except MaskedCCDBiasImageException:
            image = ccd[amp]

        ## Find spot footprints
        image -= self._bg_image(ccd, amp, *bg_reg)
        imarr = image.getImage().getArray()

        flags = afwMath.MEDIAN | afwMath.STDEVCLIP
        statistics = afwMath.makeStatistics(image, flags, ccd.stat_ctrl)
        median = statistics.getValue(afwMath.MEDIAN)
        stdev = statistics.getValue(afwMath.STDEVCLIP)

        threshold = afwDetect.Threshold(median + self.nsig*stdev)
        npix_min = 60
        if logger is not None:
            logger.info("SpotMomentFit.process_image: threshold= %s"
                        % threshold.getValue())
        fp_set = afwDetect.FootprintSet(image.getImage(), threshold, npix_min)

        ## Iterate over each footprint and calculate moments
        x0, y0 = [], [] # These are arrays to hold per-spot measurements
        prect_data = [] # This holds flattened spot footprint
        sigmax, sigmay, dn, dn_fp, chiprob = [], [], [], [], [] # Arrays to hold the moments fitting measurements
        chi2s, dofs = [], []
        maxDNs = []
        failed_curve_fits = 0
        num_fp = 0
        for fp in fp_set.getFootprints():

            num_fp +=1
            peak = [pk for pk in fp.getPeaks()][0]
            spans = fp.getSpans()
            positions = []
            zvals = []
            dn_sum = 0

            for span in spans:
                y = span.getY()
                for x in range(span.getX0(), span.getX1() + 1):
                    positions.append((x, y))
                    zvals.append(imarr[y][x])
                    dn_sum += imarr[y][x]

            try:
                # Use clipped stdev as DN error estimate for all pixels
                dn_errors = stdev*np.ones(len(positions))
                p0 = (peak.getIx(), peak.getIy(), sigma0, sigma0, dn0)
                pars, _ = scipy.optimize.leastsq(residuals, p0,
                                                     args=(positions, zvals,
                                                           dn_errors))
                sigmax.append(pars[2])
                sigmay.append(pars[3])
                dn.append(pars[4])
                x0.append(pars[0])
                y0.append(pars[1])

                dn_fp.append(dn_sum)
                chi2 = chisq(positions, zvals, x0[-1], y0[-1],
                             sigmax[-1], sigmay[-1], dn[-1], dn_errors)
                dof = fp.getNpix() - 5 # number of pars for our fit
                chiprob.append(gammaincc(dof/2., chi2/2.))
                chi2s.append(chi2)
                dofs.append(dof)
                maxDNs.append(max(zvals))

                prect_data_row = prect_values(peak,imarr, ixm=10, ixp=10, iym=10, iyp=10)
                prect_data.append(prect_data_row)
            except RuntimeError:
                failed_curve_fits +=1

        self._save_ext_data(amp, x0, y0,prect_data, sigmax, sigmay, dn, dn_fp, chiprob,
                            chi2s, dofs, maxDNs)
        self.amp_set.add(amp)

        print(num_fp)

    def numGoodFits(self, chiprob_min=0.1):
        raise NotImplementedError

    def _save_ext_data(self, amp, x0, y0, prect_data,sigmax, sigmay, dn, dn_fp, chiprob,
                        chi2s, dofs, maxDNs):
        """
        Write results from the source detection and Gaussian fitting
        to the FITS extension corresponding to the specified
        amplifier.
        """
        extname = 'Amp%02i' % amp
        try:
            #
            # Append new rows if HDU for this segment already exists.
            #
            table_hdu = self.output[extname]
            row0 = table_hdu.header['NAXIS2']
            nrows = row0 + len(x0)
            table_hdu = fitsTableFactory(table_hdu.data, nrows=nrows)
            for i in range(len(x0)):
                row = i + row0
                table_hdu.data[row]['AMPLIFIER'] = amp
                table_hdu.data[row]['XPOS'] = x0[i]
                table_hdu.data[row]['YPOS'] = y0[i]
                table_hdu.data[row]['PRECT_DATA'] = prect_data[i]
                table_hdu.data[row]['SIGMAX'] = sigmax[i]
                table_hdu.data[row]['SIGMAY'] = sigmay[i]
                table_hdu.data[row]['DN'] = dn[i]
                table_hdu.data[row]['DN_FP_SUM'] = dn_fp[i]
                table_hdu.data[row]['CHIPROB'] = chiprob[i]
                table_hdu.data[row]['CHI2'] = chi2s[i]
                table_hdu.data[row]['DOF'] = dofs[i]
                table_hdu.data[row]['MAXDN'] = maxDNs[i]
            table_hdu.name = extname
            self.output[extname] = table_hdu
        except KeyError:
            #
            # Extension for this segment does not yet exist, so add it.
            #
            colnames = ['AMPLIFIER', 'XPOS', 'YPOS', 'PRECT_DATA','SIGMAX', 'SIGMAY', 'DN',
                        'DN_FP_SUM', 'CHIPROB', 'CHI2', 'DOF', 'MAXDN']

            columns = [np.ones(len(x0))*amp, np.array(x0), np.array(y0),np.array(prect_data),
                       np.array(sigmax), np.array(sigmay),
                       np.array(dn), np.array(dn_fp),
                       np.array(chiprob), np.array(chi2s), np.array(dofs),
                       np.array(maxDNs)]
            formats = formats = ['I'] + ['E']*(2) + ['441E'] +['E']*(8)
            units = ['None', 'pixel', 'pixel', 'pixel', 'pixel',
                     'ADU', 'ADU', 'None', 'None', 'None', 'ADU',
                     'pixel', 'pixel', 'ADU']
            fits_cols = lambda coldata: [fits.Column(name=colname,
                                                     format=format,
                                                     unit=unit,
                                                     array=column)
                                         for colname, format, unit, column
                                         in coldata]
            self.output.append(fitsTableFactory(fits_cols(zip(colnames,
                                                              formats,
                                                              units,
                                                              columns))))
            self.output[-1].name = extname

    def read_spot_catalog(self, spot_catalog, chiprob_min=0.1):
        raise NotImplementedError

    def results(self, min_prob=0.1, amp=None):
        raise NotImplementedError

    def write_results(self, outfile='spot_moment_params.fits'):
        self.output[0].header['NAMPS'] = len(self.amp_set)
        fitsWriteto(self.output, outfile, overwrite=True)
