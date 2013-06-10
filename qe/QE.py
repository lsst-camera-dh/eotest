"""
@brief Class to compute medians per amplifier for a set of exposures
taken as a function of wavelength.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import glob
from collections import OrderedDict
import numpy as np
import lsst.afw.math as afwMath
import image_utils as imutils
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
    def __init__(self, results_file, pattern=None, verbose=True,
                 pd_scaling=1e-12, clobber=False):
        self.results_file = results_file
        self.verbose = verbose
        self.pd_scaling = pd_scaling
        if pattern is not None:
            self._calculate_medians(pattern, clobber)
        else:
            self._read_medians()
    def _calculate_medians(self, pattern, clobber):
        files = glob.glob(pattern)
        files.sort()
        self.medians = dict([(amp, []) for amp in imutils.allAmps])
        self.wl = []          # wavelength in nanometers
        self.exptime = []     # exposure time in seconds
        self.pd = []          # photodiode current in amps

        if os.path.isfile(self.results_file) and not clobber:
            raise RuntimeError("Results file already exists.")
        output = open(self.results_file, 'w')
        for item in files:
            if self.verbose:
                print 'processing', item
            ccd = MaskedCCD(item)
            md = Metadata(item, 1)
            exptime = md.get('EXPTIME')
            if exptime == 0:
                print "Zero exposure time in %s. Skipping." % item
                continue
            self.wl.append(md.get('MONO.WAVELENG'))
            self.exptime.append(exptime)
            if self.exptime[-1] == 0:
                raise RuntimeError("Zero exposure time in ", item)
            self.pd.append(np.abs(md.get('K_PHOT.CURRENT')*self.pd_scaling))
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
    def _read_medians(self):
        data = np.recfromtxt(self.results_file)
        data = data.transpose()
        self.wl = data[0]
        self.exptime = data[1]
        self.pd = data[2]
        self.medians = dict([(amp, col) for amp, col
                             in zip(imutils.allAmps, data[3:])])
    def calculate_QE(self, ccd_cal_file, sph_cal_file, wlscan_file,
                     gains, pixel_area=1e-10, pd_area=1e-4):
        pd_sph = PhotodiodeResponse(sph_cal_file)
        ccd_frac = CcdIllumination(wlscan_file, ccd_cal_file, sph_cal_file)
        qe = {}
        wlarrs = {}
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
        band_qe = {}
        for band in self.band_pass:
            indx = self._index(wl, self.band_pass[band])
            if len(indx[0]) != 0:
                mean_value = sum(qe[indx])/len(indx[0])
                band_qe[band] = mean_value
        return band_qe

if __name__ == '__main__':
    import pylab_plotter as plot

    ccd_cal_file = 'OD142.csv'
    sph_cal_file = 'OD143.csv'
    wlscan_file = 'WLscan.txt'

    gains = dict([(amp, 4.5) for amp in imutils.allAmps])

#    pattern = '/nfs/farm/g/lsst/u1/testData/HarvardData/112-01/final/bss50/qe/*.fits.gz'
    pattern = None
    outfile = 'med_vs_wl.txt'

    qe_data = QE_Data(outfile, pattern=pattern)
    qe_data.calculate_QE(ccd_cal_file, sph_cal_file, wlscan_file, gains)

    for i, amp in enumerate(imutils.allAmps):
        plot.curve(qe_data.wlarrs[amp], qe_data.qe[amp], oplot=i,
                   xname='wavelength (nm)', yname='QE (%)',
                   xrange=(350, 1100))
        plot.xyplot(qe_data.wlarrs[amp], qe_data.qe[amp], oplot=1)
        wl = [sum(qe_data.band_pass[band])/2. for band in qe_data.qe_band[amp]]
        plot.xyplot(wl, qe_data.qe_band[amp].values(), oplot=1, color='r')

    plot.curve(qe_data.wlarrs[1], qe_data.ccd_qe, oplot=1, color='g')
    plot.xyplot(wl, qe_data.ccd_qe_band.values(), oplot=1, color='b')
