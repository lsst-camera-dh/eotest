"""
@brief Fit 2D Gaussian to Fe55 footprints to determine the Gaussian
width of the charge dispersed signals.  For each footprint, also
compute the probability of the chi-square fit.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import numpy as np
import warnings
import itertools
import astropy.io.fits as fits
from lsst.eotest.fitsTools import fitsTableFactory, fitsWriteto
import scipy.optimize
from scipy.special import erf, gammaincc

import lsst.afw.detection as afwDetect
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath

import lsst.eotest.image_utils as imutils
from .MaskedCCD import MaskedCCD, MaskedCCDBiasImageException

import pdb

_sqrt2 = np.sqrt(2)

def psf_sigma_statistics(sigma, bins=50, range=(2, 6), frac=0.5):
    hist = np.histogram(sigma, bins=bins, range=range)
    y = hist[0]
    x = (hist[1][1:] + hist[1][:-1])/2.
    mode = x[np.where(y == max(y))[0][0]]
    my_sigma = sorted(sigma)
    median = my_sigma[int(frac*len(my_sigma))]
    mean = np.mean(sigma)
    return mode, median, mean

def pixel_integral(x, y, x0, y0, sigmax, sigmay):
    """
    Integrate 2D Gaussian centered at (x0, y0) with widths sigmax and
    sigmay over a square pixel at (x, y) with unit width.
    """
    x1, x2 = x - 0.5, x + 0.5
    y1, y2 = y - 0.5, y + 0.5

    Fx = 0.5*(erf((x2 - x0)/_sqrt2/sigmax) - erf((x1 - x0)/_sqrt2/sigmax))
    Fy = 0.5*(erf((y2 - y0)/_sqrt2/sigmay) - erf((y1 - y0)/_sqrt2/sigmay))

    return Fx*Fy

def residuals_single(pars, pos, dn, errors):
    x0, y0, sigma, DN_tot = pars
    return (dn - psf_func(pos, x0, y0, sigma, sigma, DN_tot))/errors

def residuals(pars, pos, dn, errors):
    x0, y0, sigmax, sigmay, DN_tot = pars
    return (dn - psf_func(pos, x0, y0, sigmax, sigmay, DN_tot))/errors

def psf_func(pos, x0, y0, sigmax, sigmay, DN_tot):
    """
    For a pixel location or list of pixel locations, pos, compute the
    DN per pixel for a 2D Gaussian with parameters:
    x0, y0: Gaussian mean x and y values
    sigmax, sigmay: Gaussian widths in x- and y-directions
    DN_tot: Gaussian normalization in ADU
    tie_xy: if True, then assume sigmax=sigmay in the fit.
    """
    return DN_tot*np.array([pixel_integral(x[0], x[1], x0, y0,
                                           sigmax, sigmay) for x in pos])

def chisq(pos, dn, x0, y0, sigmax, sigmay, dn_fit, dn_errors):
    "The chi-square of the fit of the data to psf_func."
    return sum((psf_func(pos, x0, y0, sigmax, sigmay, dn_fit)
                - np.array(dn))**2/dn_errors**2)

def p9_values(peak, imarr, x0, y0, sigmax, sigmay, DN_tot):
    x5, y5 = peak.getIx(), peak.getIy()
    pos = [(x5 + dyx[1], y5 + dyx[0]) for dyx in
           itertools.product((-1, 0, 1), (-1, 0, 1))]
    p9_data = np.array([imarr[y][x] for x, y in pos])
    p9_model = psf_func(pos, x0, y0, sigmax, sigmay, DN_tot)
    return p9_data, p9_model

def prect_values(peak, imarr, ixm=3, ixp=21, iym=3, iyp=3):
    xpeak, ypeak = peak.getIx(), peak.getIy()
    # nb. imarr includes overscan
    yimsiz,ximsiz = imarr.shape
    # if we are too close to an edge, just return all 0's
    if ypeak-iym<0 or xpeak-ixm<0 or ypeak+iyp+1>yimsiz or xpeak+ixp+1>ximsiz:
        prect_data = np.zeros([iyp+iym+1, ixp+ixm+1])
    else:
        # store a region around peak pixel [-3,21] in x and [-3,3] in y
        prect_data = imarr[ypeak-iym:ypeak+iyp+1, xpeak-ixm:xpeak+ixp+1]

    # data is stored as a normal 2-d array and then flattened
    prect_data_flat = prect_data.flatten()
    if prect_data_flat.shape[0] == 0:
        pdb.set_trace()
    return prect_data_flat

class PsfGaussFit(object):
    def __init__(self, nsig=3, min_npix=None, max_npix=20, gain_est=2,
                 fit_xy=True, outfile=None):
        """
        nsig is the threshold in number of clipped stdev above median.

        If fit_xy == True, then the Gaussian widths in the x- and
        y-directions are fit separately, otherwise a single Gaussian
        width is assumed for both directions.

        min_npix is the minimum number of pixels to be used in the
        4- or 5-parameter fit.  If min_npix is None, then it is set to
        5 for a 4-parameter fit (fit_xy==False), otherwise it is set
        to 6 for a 5-parameter fit (fit_xy==True).
        """
        self.nsig = nsig
        if fit_xy:
            self.npars = 5
        else:
            self.npars = 4
        if min_npix is None:
            self.min_npix = self.npars + 1
        else:
            self.min_npix = min_npix
        self.max_npix = max_npix
        self.sigmax, self.sigmay = [], []
        self.dn, self.dn_fp, self.chiprob = [], [], []
        self.amp = []
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
        """
        Process a segment and accumulate the fit results for each
        charge cluster.  The dn0 and sigma0 parameters are the
        starting values used for each fit.
        """
        try:
            image = ccd.bias_subtracted_image(amp, fit_order=oscan_fit_order)
        except MaskedCCDBiasImageException:
            print "DM stack error encountered when generating bias image "
            print "from inferred overscan region."
            print "Skipping bias subtraction."
            image = ccd[amp]

        image -= self._bg_image(ccd, amp, *bg_reg)
        imarr = image.getImage().getArray()

        flags = afwMath.MEDIAN | afwMath.STDEVCLIP
        statistics = afwMath.makeStatistics(image, flags, ccd.stat_ctrl)
        median = statistics.getValue(afwMath.MEDIAN)
        stdev = statistics.getValue(afwMath.STDEVCLIP)

        threshold = afwDetect.Threshold(median + self.nsig*stdev)
        if logger is not None:
            logger.info("PsfGaussFit.process_image: threshold= %s"
                        % threshold.getValue())
        fpset = afwDetect.FootprintSet(image, threshold)

        x0, y0 = [], []
        sigmax, sigmay, dn, dn_fp, chiprob = [], [], [], [], []
        chi2s, dofs = [], []
        maxDNs = []
        xpeak, ypeak = [], []
        p9_data = []
        p9_model = []
        prect_data = []
        failed_curve_fits = 0
        num_fp = 0
        for fp in fpset.getFootprints():
            if fp.getNpix() < self.min_npix or fp.getNpix() > self.max_npix:
                continue
            num_fp += 1
            spans = fp.getSpans()
            positions = []
            zvals = []
            peak = [pk for pk in fp.getPeaks()][0]
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
                if self.npars == 5:
                    p0 = (peak.getIx(), peak.getIy(), sigma0, sigma0, dn0)
                    pars, _ = scipy.optimize.leastsq(residuals, p0,
                                                     args=(positions, zvals,
                                                           dn_errors))
                    sigmax.append(pars[2])
                    sigmay.append(pars[3])
                    dn.append(pars[4])
                else:
                    p0 = (peak.getIx(), peak.getIy(), sigma0, dn0)
                    pars, _ = scipy.optimize.leastsq(residuals_single, p0,
                                                     args=(positions, zvals,
                                                           dn_errors))
                    sigmax.append(pars[2])
                    sigmay.append(pars[2])
                    dn.append(pars[3])
                x0.append(pars[0])
                y0.append(pars[1])
                dn_fp.append(dn_sum)
                chi2 = chisq(positions, zvals, x0[-1], y0[-1],
                             sigmax[-1], sigmay[-1], dn[-1], dn_errors)
                dof = fp.getNpix() - self.npars
                chiprob.append(gammaincc(dof/2., chi2/2.))
                chi2s.append(chi2)
                dofs.append(dof)
                maxDNs.append(max(zvals))
                try:
                    p9_data_row, p9_model_row \
                        = p9_values(peak, imarr, x0[-1], y0[-1], sigmax[-1],
                                    sigmay[-1], dn[-1])
                    prect_data_row = prect_values(peak,imarr)
                    p9_data.append(p9_data_row)
                    p9_model.append(p9_model_row)
                    prect_data.append(prect_data_row)
                    xpeak.append(peak.getIx())
                    ypeak.append(peak.getIy())
                except IndexError:
                    [item.pop() for item in (x0, y0, sigmax, sigmay,
                                             dn, dn_fp, chiprob,
                                             chi2s, dofs, maxDNs)]
            except RuntimeError:
                failed_curve_fits += 1
                pass
        if logger is not None:
            logger.info("Number of footprints fitted: %i" % num_fp)
            if failed_curve_fits > 0:
                logger.info("Failed scipy.optimize.leastsq calls: %s"
                            % failed_curve_fits)
        self._save_ext_data(amp, x0, y0, sigmax, sigmay, dn, dn_fp, chiprob,
                            chi2s, dofs, maxDNs, xpeak, ypeak,
                            np.array(p9_data), np.array(p9_model),
                            np.array(prect_data))
        self.amp_set.add(amp)
        self.sigmax.extend(sigmax)
        self.sigmay.extend(sigmay)
        self.dn.extend(dn)
        self.dn_fp.extend(dn_fp)
        self.chiprob.extend(chiprob)
        self.amp.extend(np.ones(len(sigmax))*amp)
    def numGoodFits(self, chiprob_min=0.1):
        chiprob = np.array(self.chiprob)
        amps = np.sort(np.unique(np.array(self.amp)))
        my_numGoodFits = dict([(amp, 0) for amp in imutils.allAmps()])
        for amp in amps:
            indx = np.where((self.chiprob > chiprob_min) & (self.amp == amp))
            my_numGoodFits[amp] = len(indx[0])
        return my_numGoodFits
    def _save_ext_data(self, amp, x0, y0, sigmax, sigmay, dn, dn_fp, chiprob,
                       chi2s, dofs, maxDNs, xpeak, ypeak, p9_data, p9_model, prect_data):
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
                table_hdu.data[row]['SIGMAX'] = sigmax[i]
                table_hdu.data[row]['SIGMAY'] = sigmay[i]
                table_hdu.data[row]['DN'] = dn[i]
                table_hdu.data[row]['DN_FP_SUM'] = dn_fp[i]
                table_hdu.data[row]['CHIPROB'] = chiprob[i]
                table_hdu.data[row]['CHI2'] = chi2s[i]
                table_hdu.data[row]['DOF'] = dofs[i]
                table_hdu.data[row]['MAXDN'] = maxDNs[i]
                table_hdu.data[row]['XPEAK'] = xpeak[i]
                table_hdu.data[row]['YPEAK'] = ypeak[i]
                table_hdu.data[row]['P9_DATA'] = p9_data[i]
                table_hdu.data[row]['P9_MODEL'] = p9_model[i]
                table_hdu.data[row]['PRECT_DATA'] = prect_data[i]
            table_hdu.name = extname
            self.output[extname] = table_hdu
        except KeyError:
            #
            # Extension for this segment does not yet exist, so add it.
            #
            colnames = ['AMPLIFIER', 'XPOS', 'YPOS', 'SIGMAX', 'SIGMAY', 'DN',
                        'DN_FP_SUM', 'CHIPROB', 'CHI2', 'DOF', 'MAXDN',
                        'XPEAK', 'YPEAK', 'P9_DATA', 'P9_MODEL', 'PRECT_DATA']
            columns = [np.ones(len(x0))*amp, np.array(x0), np.array(y0),
                       np.array(sigmax), np.array(sigmay),
                       np.array(dn), np.array(dn_fp),
                       np.array(chiprob), np.array(chi2s), np.array(dofs),
                       np.array(maxDNs), np.array(xpeak), np.array(ypeak),
                       np.array(p9_data), np.array(p9_model), np.array(prect_data)]
            formats = ['I'] + ['E']*(len(columns)-6) + ['I']*2 + ['9E']*2 + ['175E']
            units = ['None', 'pixel', 'pixel', 'pixel', 'pixel',
                     'ADU', 'ADU', 'None', 'None', 'None', 'ADU',
                     'pixel', 'pixel', 'ADU', 'ADU', 'ADU']
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
    def read_fe55_catalog(self, psf_catalog, chiprob_min=0.1):
        catalog = fits.open(psf_catalog)
        for attr in 'sigmax sigmay dn dn_fp_sum chiprob amp'.split():
            exec('self.%s = np.array((), dtype=float)' % attr)
        for amp in imutils.allAmps(psf_catalog):
            extname = 'Amp%02i' % amp
            chiprob = catalog[extname].data.field('CHIPROB')
            index = np.where(chiprob > chiprob_min)
            self.chiprob = np.concatenate((self.chiprob, chiprob[index]))
            self.amp = np.concatenate((self.amp, np.ones(len(index[0]))*amp))
            for attr in 'sigmax sigmay dn dn_fp_sum'.split():
                command = 'self.%(attr)s = np.concatenate((self.%(attr)s, catalog["%(extname)s"].data.field("%(attr)s")[index]))' % locals()
                exec(command)
            self.dn_fp = self.dn_fp_sum
    def results(self, min_prob=0.1, amp=None):
        """
        Return sigmax, sigmay, dn, chiprob for chiprob > min_prob for
        specified amp.
        """
        chiprob = np.array(self.chiprob, dtype=np.float)
        amps = np.array(self.amp, dtype=np.int)
        if amp is not None:
            indx = np.where((chiprob > min_prob) & (amps == amp))
        else:
            indx = np.where(chiprob > min_prob)
        my_results = {}
        my_results['sigmax'] = np.array(self.sigmax, dtype=np.float)[indx]
        my_results['sigmay'] = np.array(self.sigmay, dtype=np.float)[indx]
        my_results['dn'] = np.array(self.dn, dtype=np.float)[indx]
        my_results['dn_fp'] = np.array(self.dn_fp, dtype=np.float)[indx]
        my_results['chiprob'] = chiprob[indx]
        my_results['amps'] = amps[indx]
        return my_results
    def write_results(self, outfile='fe55_psf_params.fits'):
        self.output[0].header['NAMPS'] = len(self.amp_set)
        fitsWriteto(self.output, outfile, clobber=True, checksum=True)

if __name__ == '__main__':
    import os
    import pylab_plotter as plot
    plot.pylab.ion()

    infile = os.path.join(os.environ['EOTESTDIR'],
                          'work/sensorData/000-00/fe55/debug/000-00_fe55_fe55_00_debug.fits')
    #infile = 'fe55_0060s_000.fits'
    outfile = '000-00_fe55_psf.fits'
    nsig = 2
    #fit_xy = False
    fit_xy = True

    ccd = MaskedCCD(infile)

    fitter = PsfGaussFit(nsig=nsig, fit_xy=fit_xy)
    for amp in ccd.keys()[:2]:
        print 'processing amp:', amp
        fitter.process_image(ccd, amp)
    fitter.write_results(outfile)

    fitter = PsfGaussFit(nsig=nsig, outfile=outfile, fit_xy=fit_xy)
    for amp in ccd.keys()[2:]:
        print "processing amp:", amp
        fitter.process_image(ccd, amp)
    fitter.write_results(outfile)

    results = fitter.results()

    flags = afwMath.MEDIAN | afwMath.STDEVCLIP

    stats = afwMath.makeStatistics(results['sigmax'], flags)
    median = stats.getValue(afwMath.MEDIAN)
    stdev = stats.getValue(afwMath.STDEVCLIP)
    plot.histogram(results['sigmax'], xname='Fitted sigma values',
                   xrange=(median-3*stdev, median+3*stdev))
    plot.histogram(results['sigmay'], oplot=1, color='r',
                   xrange=(median-3*stdev, median+3*stdev))

    plot.histogram(results['dn'], xname='Fitted DN values', xrange=(250, 450))

    plot.xyplot(results['chiprob'], results['sigmax'],
                xname='chi-square prob.', yname='sigma', ylog=1)
                
    plot.xyplot(results['chiprob'], results['sigmay'], oplot=1, color='r')

    plot.xyplot(results['dn'], results['dn_fp'], xname='DN (fitted value)',
                yname='DN (footprint sum)', xrange=(0, max(results['dn_fp'])),
                yrange=(0, max(results['dn_fp'])))
