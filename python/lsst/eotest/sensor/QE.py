"""
@brief Class to compute medians per amplifier for a set of exposures
taken as a function of wavelength.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import sys
import glob
from collections import OrderedDict
import numpy as np
import pyfits
from lsst.eotest.pyfitsTools import pyfitsTableFactory, pyfitsWriteto
import lsst.eotest.image_utils as imutils
import pylab_plotter as plot
from MaskedCCD import MaskedCCD
from PhotodiodeResponse import PhotodiodeResponse, CcdIllumination, \
     Interpolator

import lsst.afw.image as afwImage
import lsst.afw.math as afwMath

planck = 6.626e-34      # Planck constant in SI
clight = 2.99792458e8   # speed of light in m/s

class QE_Data(object):
    """Class to contain sequences of image medians (in DN) indexed by
    amplifier for a set of exposures made at a series of incident
    illumination wavelengths.  The wavelength, exposure time, and
    photodiode current are also recorded for each exposure."""
    band_pass = OrderedDict()
    band_pass['u'] = (321, 391)
    band_pass['g'] = (402, 552)
    band_pass['r'] = (552, 691)
    band_pass['i'] = (691, 818)
    band_pass['z'] = (818, 922)
    band_pass['y'] = (930, 1070)
    band_wls = np.array([sum(band_pass[b])/2. for b in band_pass.keys()])
    def __init__(self, verbose=True, pd_scaling=1e-9, logger=None):
        self.verbose = verbose
        self.pd_scaling = pd_scaling
        self.logger = logger
    def read_medians(self, medians_file):
        data = np.recfromtxt(medians_file)
        data = data.transpose()
        self.wl = data[0]
        self.exptime = data[1]
        self.pd = data[2]
        self.medians = dict([(amp, col) for amp, col
                             in zip(imutils.allAmps, data[3:])])
    def calculate_medians(self, infiles, outfile, mask_files=(),
                          clobber=False):
        files = sorted([x for x in infiles])
        self.medians = dict([(amp, []) for amp in imutils.allAmps])
        self.wl = []          # wavelength in nanometers
        self.exptime = []     # exposure time in seconds
        self.pd = []          # photodiode current in amps

        if os.path.isfile(outfile) and not clobber:
            raise RuntimeError("Output file for image medians already exists.")
        output = open(outfile, 'w')
        for item in files:
            if self.verbose:
                if self.logger is not None:
                    self.logger.info('processing %s' % item)
                else:
                    print 'processing', item
            ccd = MaskedCCD(item, mask_files=mask_files)
            md = imutils.Metadata(item, 1)
            wl = md.get('MONOWL')
            exptime = md.get('EXPTIME')
            if exptime == 0:
                line = "Zero exposure time in %s. Skipping." % item
                if self.logger is not None:
                    self.logger.info(line + '\n')
                else:
                    print line 
                continue
            self.wl.append(wl)
            self.exptime.append(exptime)
            if self.exptime[-1] == 0:
                raise RuntimeError("Zero exposure time in ", item)
            try:
                self.pd.append(np.abs(md.get('MONDIODE')*self.pd_scaling))
            except TypeError, e:
                if self.logger is not None:
                    self.logger.info('%s\n' % e)
                continue
            output.write('%.3f' % self.wl[-1])
            output.write('  %.3e' % self.exptime[-1])
            output.write('  %.3e' % self.pd[-1])
            for amp in imutils.allAmps:
                im = ccd.unbiased_and_trimmed_image(amp)
                value = afwMath.makeStatistics(im, afwMath.MEDIAN,
                                               ccd.stat_ctrl).getValue()
                self.medians[amp].append(value)
                output.write('  %.5f' % self.medians[amp][-1])
            output.write('\n')
            output.flush()
        output.close()
        self.wl = np.array(self.wl)
        self.exptime = np.array(self.exptime)
        self.pd = np.abs(np.array(self.pd))
    def incidentPower_e2v(self):
        # For e2v data MONDIODE = incident power/area = nW/cm**2 so
        # just need to multiply by times nominal pixel area in
        # cm**2.
        power = []
        pixel_area = 1e-3*1e-3 # (10 microns)**2
        for pd_current in self.pd:
            power.append(pd_current*pixel_area)
        self.power = np.array(power, dtype=np.float)
    def incidentPower(self, pd_ratio_file, pixel_area=1e-10):
        # Calibration diode collecting area.
        pd_area = float(open(pd_ratio_file).readline().split('=')[1])
        if self.logger is not None:
            self.logger.info("Using pd_ratio_file: %s" % pd_ratio_file)
            self.logger.info("pd_area = %.2e mm^2" % pd_area)
        # Incident power per pixel (J/s)
        data = np.recfromtxt(pd_ratio_file, skip_header=2,
                             names='monowl, sens, ccdfrac')
        sensitivity = Interpolator(data['monowl'], data['sens'])
        ccd_frac = Interpolator(data['monowl'], data['ccdfrac'])
        power = []
        for pd_current, wl_nm in zip(self.pd, self.wl):
            power.append(pd_current/(ccd_frac(wl_nm)*sensitivity(wl_nm))
                         *pixel_area/pd_area)
        self.power = np.array(power, dtype=np.float)
    def calculate_QE(self, gains):
        qe = OrderedDict()
        wlarrs = OrderedDict()
        for amp in imutils.allAmps:
            qe[amp] = []
            wlarrs[amp] = []
            for i, wl_nm in enumerate(self.wl):
                try:
                    wl = wl_nm*1e-9
                    hnu = planck*clight/wl   # photon energy (J)
                    #
                    # Incident illumination power per pixel (J/s)
                    #
                    power = self.power[i]
                    #
                    # Median number of photo-electrons per pixel
                    #
                    Ne = self.medians[amp][i]*gains[amp]
                    #
                    # Number of incident photons per pixel.
                    #
                    nphot = self.exptime[i]*power/hnu
                    qe_value = Ne/nphot
                    if qe_value > 10.:
                        # Skip obviously out-of-range values
                        raise RuntimeError("QE value > 10")
                    qe[amp].append(100*qe_value)
                    wlarrs[amp].append(wl_nm)
                except RuntimeError:
                    # Expect a RuntimeError if the requested wavelength is
                    # outside of the photo-diode calibration or outside of
                    # the wavelength scan that compares the illumination at
                    # the CCD vs at the integrating sphere.  Skip these
                    # cases.
                    continue
            qe[amp] = np.abs(np.array(qe[amp]))
            wlarrs[amp] = np.array(wlarrs[amp])
        self.qe = qe
        self.wlarrs = wlarrs
        self.qe_band = dict([(amp, {}) for amp in self.qe])
        for amp in self.qe:
            self.qe_band[amp] = self.compute_band_qe(self.wlarrs[amp],
                                                     self.qe[amp])
        self.ccd_qe = sum(self.qe.values())/len(self.qe.values())
        self.ccd_qe_band = self.compute_band_qe(self.wlarrs[1], self.ccd_qe)
    def _index(self, wl, band_pass):
        indx = np.where((wl >= band_pass[0]) & (wl <= band_pass[1]))
        return indx
    def compute_band_qe(self, wl, qe):
        band_qe = OrderedDict()
        for band in self.band_pass:
            indx = self._index(wl, self.band_pass[band])
            if len(indx[0]) != 0:
                mean_value = sum(qe[indx])/len(indx[0])
                band_qe[band] = mean_value
        return band_qe
    def write_fits_tables(self, outfile, clobber=True):
        amps = self.qe.keys()
        colnames = ['WAVELENGTH']
        colnames.extend(['AMP%02i' % i for i in amps])
        colnames.append('DEVICE_MEAN')

        columns = [self.wlarrs[1]]
        columns.extend([self.qe[i] for i in amps])
        columns.append(self.ccd_qe)

        formats = ["E"]*len(columns)
        units = ["nm"]
        units.extend(["e-/photon %"]*(len(columns)-1))

        fits_cols = lambda coldata : [pyfits.Column(name=colname,
                                                    format=format,
                                                    unit=unit, array=column) 
                                      for colname, format, unit, column
                                      in coldata]

        HDUList = pyfits.HDUList()
        HDUList.append(pyfits.PrimaryHDU())
        HDUList.append(pyfitsTableFactory(fits_cols(zip(colnames, formats,
                                                        units, columns))))
        HDUList[-1].name = 'QE_CURVES'

        columns = [self.ccd_qe_band.keys()]
        columns.extend([np.array(self.qe_band[amp].values()) for amp in amps])
        columns.append(np.array(self.ccd_qe_band.values()))

        colnames[0] = 'BAND'
        formats[0] = '2A'
        units[0] = None

        HDUList.append(pyfitsTableFactory(fits_cols(zip(colnames, formats,
                                                        units, columns))))
        HDUList[-1].name = 'QE_BANDS'
        pyfitsWriteto(HDUList, outfile, clobber=clobber)
    def plot_curves(self, outfile=None, interactive=False):
        if interactive:
            plot.pylab.ion()
        else:
            plot.pylab.ioff()
        indx = np.argsort(self.wlarrs[1])
        for i, amp in enumerate(imutils.allAmps):
            plot.curve(self.wlarrs[amp][indx], self.qe[amp][indx],
                       oplot=i, xname='wavelength (nm)', yname='QE (%)',
                       xrange=(350, 1100))
            plot.xyplot(self.wlarrs[amp][indx], self.qe[amp][indx], oplot=1)
            wl = [sum(self.band_pass[band])/2. for band in self.qe_band[amp]]
            plot.xyplot(wl, self.qe_band[amp].values(), oplot=1, color='r')
        plot.curve(self.wlarrs[1][indx], self.ccd_qe[indx], oplot=1,
                   color='g')
        plot.xyplot(wl, self.ccd_qe_band.values(), oplot=1, color='b')
        if outfile is not None:
            plot.save(outfile)

if __name__ == '__main__':
    ccd_cal_file = 'OD142.csv'
    sph_cal_file = 'OD143.csv'
    wlscan_file = 'WLscan.txt'

    gains = dict([(amp, 2.5) for amp in imutils.allAmps])
    infiles = glob.glob('/nfs/farm/g/lsst/u1/testData/HarvardData/112-01/final/bss70/qe/112_01_qe_[0-9]*.fits.gz')
#    medians_file = 'med_vs_wl.txt'
    medians_file = 'med_vs_wl_bss70.txt'

    qe_data = QE_Data()
#    qe_data.calculate_medians(infiles, medians_file)
    qe_data.read_medians(medians_file)
    qe_data.incidentPower(ccd_cal_file, sph_cal_file, wlscan_file)
    qe_data.calculate_QE(gains)
    qe_data.write_fits_tables('qe_tables.fits')
    qe_data.plot_curves(interactive=True)
