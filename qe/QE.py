"""
@brief Class to compute medians per amplifier for a set of exposures
taken as a function of wavelength.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import glob
from collections import OrderedDict
import numpy as np
import pyfits
import lsst.afw.math as afwMath
import image_utils as imutils
import pylab_plotter as plot
from MaskedCCD import MaskedCCD, Metadata
from PhotodiodeResponse import PhotodiodeResponse, CcdIllumination

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
    def __init__(self, verbose=True, pd_scaling=1e-9):
        self.verbose = verbose
        self.pd_scaling = pd_scaling
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
                print 'processing', item
            ccd = MaskedCCD(item, mask_files=mask_files)
            md = Metadata(item, 1)
            exptime = md.get('EXPTIME')
            if exptime == 0:
                print "Zero exposure time in %s. Skipping." % item
                continue
            self.wl.append(md.get('MONOWL'))
            self.exptime.append(exptime)
            if self.exptime[-1] == 0:
                raise RuntimeError("Zero exposure time in ", item)
            self.pd.append(np.abs(md.get('MONDIODE')*self.pd_scaling))
            output.write('%.3f' % self.wl[-1])
            output.write('  %.3e' % self.exptime[-1])
            output.write('  %.3e' % self.pd[-1])
            for amp in imutils.allAmps:
                im = imutils.unbias_and_trim(ccd[amp])
                value = afwMath.makeStatistics(im, afwMath.MEDIAN).getValue()
                self.medians[amp].append(value)
                output.write('  %.5f' % self.medians[amp][-1])
            output.write('\n')
            output.flush()
        output.close()
        self.wl = np.array(self.wl)
        self.exptime = np.array(self.exptime)
        self.pd = np.abs(np.array(self.pd))
    def calculate_QE(self, ccd_cal_file, sph_cal_file, wlscan_file,
                     gains, pixel_area=1e-10, pd_area=1e-4):
        pd_sph = PhotodiodeResponse(sph_cal_file)
        ccd_frac = CcdIllumination(wlscan_file, ccd_cal_file, sph_cal_file)
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
                    power = ((self.pd[i]/pd_sph(wl_nm))*ccd_frac(wl_nm)
                             *(pixel_area/pd_area))
                    #
                    # Median number of photo-electrons per pixel
                    #
                    Ne = self.medians[amp][i]*gains[amp]
                    #
                    # Number of incident photons per pixel.
                    #
                    nphot = self.exptime[i]*power/hnu
                    qe_value = Ne/nphot
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
        HDUList.append(pyfits.new_table(fits_cols(zip(colnames, formats,
                                                      units, columns))))
        HDUList[-1].name = 'QE_CURVES'

        columns = [self.ccd_qe_band.keys()]
        columns.extend([np.array(self.qe_band[amp].values()) for amp in amps])
        columns.append(np.array(self.ccd_qe_band.values()))

        colnames[0] = 'BAND'
        formats[0] = '2A'
        units[0] = None

        HDUList.append(pyfits.new_table(fits_cols(zip(colnames, formats,
                                                      units, columns))))
        HDUList[-1].name = 'QE_BANDS'
        HDUList.writeto(outfile, clobber=clobber)
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
    qe_data.calculate_QE(ccd_cal_file, sph_cal_file, wlscan_file, gains)
    qe_data.write_fits_tables('qe_tables.fits')
    qe_data.plot_curves(interactive=True)
