"""
@brief Fit 2D Gaussian to Fe55 footprints to determine the Gaussian
width of the charge dispersed signals.  For each footprint, also
compute the probability of the chi-square fit.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import numpy as np
import pyfits
import scipy.optimize
from scipy.special import erf, gammaincc

import lsst.afw.detection as afwDetect
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath

import image_utils as imutils
from MaskedCCD import MaskedCCD

_sqrt2 = np.sqrt(2)

def pixel_integral(x, y, x0, y0, sigma):
    """
    Integrate 2D Gaussian centered at (x0, y0) with width sigma over a
    square pixel at (x, y) with unit width.
    """
    x1, x2 = x - 0.5, x + 0.5
    y1, y2 = y - 0.5, y + 0.5

    Fx = 0.5*(erf((x2 - x0)/_sqrt2/sigma) - erf((x1 - x0)/_sqrt2/sigma))
    Fy = 0.5*(erf((y2 - y0)/_sqrt2/sigma) - erf((y1 - y0)/_sqrt2/sigma))
    
    return Fx*Fy

def psf_func(pos, *args):
    """
    For a pixel location or list of pixel locations, pos, compute the
    DN for a 2D Gaussian with parameters args.
    """
    x0, y0, sigma, DN = args
    if type(pos) == type([]):
        return DN*np.array([pixel_integral(x[0], x[1], x0, y0, sigma) 
                            for x in pos])
    return DN*pixel_integral(x[0], x[1], x0, y0, sigma)

def chisq(pos, dn, args, dn_errors):
    "The chi-square of the fit of the data to psf_func."
    return sum((psf_func(pos, *tuple(args)) - np.array(dn))**2/dn_errors**2)

class PsfGaussFit(object):
    def __init__(self, nsig=3, min_npix=5, gain_est=2, outfile=None):
        """
        nsig is the threshold in number of clipped stdev above median.
        min_npix is the minimum number of pixels to be used in the
        4-parameter fit.
        """
        self.nsig = nsig
        self.min_npix = min_npix
        self.sigma, self.dn, self.dn_fp, self.chiprob = [], [], [], []
        self.amp = []
        self.outfile = outfile
        if outfile is None:
            self.output = pyfits.HDUList()
            self.output.append(pyfits.PrimaryHDU())
        else:
            # Append new data to existing file.
            self.output = pyfits.open(self.outfile)
    def _bg_image(self, ccd, amp, nx, ny):
        "Compute background image based on clipped local mean."
        bg_ctrl = afwMath.BackgroundControl(nx, ny, ccd.stat_ctrl)
        bg = afwMath.makeBackground(ccd[amp], bg_ctrl)
        return bg.getImageF()
    def process_image(self, ccd, amp, sigma0=0.36, dn0=1590./5.,
                      bg_reg=(10, 10)):
        """
        Process a segment and accumulate the fit results for each
        charge cluster.  The dn0 and sigma0 parameters are the
        starting values used for each fit.
        """
        image = ccd.bias_subtracted_image(amp)
        image -= self._bg_image(ccd, amp, *bg_reg)
        imarr = image.getImage().getArray()

        flags = afwMath.MEDIAN | afwMath.STDEVCLIP
        statistics = afwMath.makeStatistics(image, flags, ccd.stat_ctrl)
        median = statistics.getValue(afwMath.MEDIAN)
        stdev = statistics.getValue(afwMath.STDEVCLIP)

        threshold = afwDetect.Threshold(median + self.nsig*stdev)
        fpset = afwDetect.FootprintSet(image, threshold)

        x0, y0 = [], []
        sigma, dn, dn_fp, chiprob = [], [], [], []
        chi2s, dofs = [], []
        maxDNs = []
        for fp in fpset.getFootprints():
            if fp.getNpix() < self.min_npix:
                continue
            spans = fp.getSpans()
            positions = []
            zvals = []
            peak = [pk for pk in fp.getPeaks()][0]
            p0 = (peak.getIx(), peak.getIy(), sigma0, dn0)
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
                pars, _ = scipy.optimize.curve_fit(psf_func, positions, 
                                                   zvals, p0=p0,
                                                   sigma=dn_errors)
                x0.append(pars[0])
                y0.append(pars[1])
                sigma.append(pars[2])
                dn.append(pars[3])
                dn_fp.append(dn_sum)
                chi2 = chisq(positions, zvals, pars, dn_errors)
                dof = fp.getNpix() - 4
                chiprob.append(gammaincc(dof/2., chi2/2.))
                chi2s.append(chi2)
                dofs.append(dof)
                maxDNs.append(max(zvals))
            except RuntimeError:
                pass
        self._save_ext_data(amp, x0, y0, sigma, dn, dn_fp, chiprob,
                            chi2s, dofs, maxDNs)
        self.sigma.extend(sigma)
        self.dn.extend(dn)
        self.dn_fp.extend(dn_fp)
        self.chiprob.extend(chiprob)
        self.amp.extend(np.ones(len(sigma))*amp)
    def _save_ext_data(self, amp, x0, y0, sigma, dn, dn_fp, chiprob,
                       chi2s, dofs, maxDNs):
        """
        Write results from the source detection and Gaussian fitting
        to the FITS extension corresponding to the specified
        amplifier.
        """
        extname = 'Segment%s' % imutils.channelIds[amp]
        try:
            #
            # Append new rows if HDU for this segment already exists.
            #
            table_hdu = self.output[extname]
            row0 = table_hdu.header['NAXIS2']
            nrows = row0 + len(x0)
            table_hdu = pyfits.new_table(table_hdu.data, nrows=nrows)
            for i in range(len(x0)):
                row = i + row0
                table_hdu.data[row]['AMPLIFIER'] = amp
                table_hdu.data[row]['XPOS'] = x0[i]
                table_hdu.data[row]['YPOS'] = y0[i]
                table_hdu.data[row]['SIGMA'] = sigma[i]
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
            colnames = ['AMPLIFIER', 'XPOS', 'YPOS', 'SIGMA', 'DN',
                        'DN_FP_SUM', 'CHIPROB', 'CHI2', 'DOF', 'MAXDN']
            columns = [np.ones(len(x0))*amp, np.array(x0), np.array(y0),
                       np.array(sigma), np.array(dn), np.array(dn_fp),
                       np.array(chiprob), np.array(chi2s), np.array(dofs),
                       np.array(maxDNs)]
            formats = ['I'] + ['E']*(len(columns)-1)
            units = ['None', 'pixel', 'pixel', 'pixel', 'ADU', 'ADU', 'None',
                     'None', 'None', 'ADU']
            fits_cols = lambda coldata : [pyfits.Column(name=colname,
                                                        format=format,
                                                        unit=unit,
                                                        array=column)
                                          for colname, format, unit, column
                                          in coldata]
            self.output.append(pyfits.new_table(fits_cols(zip(colnames,
                                                              formats,
                                                              units,
                                                              columns))))
            self.output[-1].name = extname
    def results(self, min_prob=0.1):
        """
        Return sigma, dn, chiprob for chiprob > min_prob.
        """
        sigma = np.array(self.sigma, dtype=np.float)
        dn = np.array(self.dn, dtype=np.float)
        dn_fp = np.array(self.dn_fp, dtype=np.float)
        chiprob = np.array(self.chiprob, dtype=np.float)
        amp = np.array(self.amp, dtype=np.int)
        indx = np.where(chiprob > min_prob)
        return sigma[indx], dn[indx], dn_fp[indx], chiprob[indx], amp[indx]
    def write_results(self, outfile='fe55_psf_params.fits'):
        self.output.writeto(outfile, clobber=True, checksum=True)

if __name__ == '__main__':
    import pylab_plotter as plot
    plot.pylab.ion()

    infile = 'work/sensorData/000-00/fe55/debug/000-00_fe55_fe55_00_debug.fits'
    #infile = 'fe55_0060s_000.fits'
    outfile = '000-00_fe55_psf.fits'
    nsig = 2

    ccd = MaskedCCD(infile)

    fitter = PsfGaussFit(nsig=nsig)
    for amp in imutils.allAmps[:2]:
        print 'processing amp:', amp
        fitter.process_image(ccd, amp)
    fitter.write_results(outfile)

    fitter = PsfGaussFit(nsig=nsig, outfile=outfile)
    for amp in imutils.allAmps[2:]:
        print "processing amp:", amp
        fitter.process_image(ccd, amp)
    fitter.write_results(outfile)

    sigma, dn, dn_fp, chiprob, amp = fitter.results()

    flags = afwMath.MEDIAN | afwMath.STDEVCLIP

    stats = afwMath.makeStatistics(sigma, flags)
    median = stats.getValue(afwMath.MEDIAN)
    stdev = stats.getValue(afwMath.STDEVCLIP)
    plot.histogram(sigma, xname='Fitted sigma values',
                   xrange=(median-3*stdev, median+3*stdev))

    plot.histogram(dn, xname='Fitted DN values', xrange=(250, 450))

    plot.xyplot(chiprob, sigma, xname='chi-square prob.', yname='sigma')

    plot.xyplot(dn, dn_fp, xname='DN (fitted value)',
                yname='DN (footprint sum)')
