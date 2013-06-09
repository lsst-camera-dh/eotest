"""
@brief Class to compute medians per amplifier for a set of exposures
taken as a function of wavelength.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import glob
import numpy as np
import lsst.afw.math as afwMath
import image_utils as imutils
from MaskedCCD import MaskedCCD, Metadata
from PhotodiodeResponse import PhotodiodeResponse, CcdIllumination

class QE_Data(object):
    """Class to contain sequences of image medians (in DN) indexed by
    amplifier for a set of exposures made at a series of incident
    illumination wavelengths.  The wavelength, exposure time, and
    photodiode current are also recorded for each exposure."""
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

if __name__ == '__main__':
    import pylab_plotter as plot

    planck = 6.626e-34      # Planck constant in SI
    clight = 2.99792458e8   # speed of light in m/s

    ccd_cal_file = 'OD142.csv'
    sph_cal_file = 'OD143.csv'
    wlscan_file = 'WLscan.txt'

    gains = dict([(amp, 5) for amp in imutils.allAmps])

#    pattern = '/nfs/farm/g/lsst/u1/testData/HarvardData/112-01/final/bss50/qe/*.fits.gz'
    pattern = None
    outfile = 'med_vs_wl.txt'
    qe_data = QE_Data(outfile, pattern=pattern)

    pd_sph = PhotodiodeResponse(sph_cal_file)
    ccd_frac = CcdIllumination(wlscan_file, ccd_cal_file, sph_cal_file)

    pixel_area = (1e-5)**2         # pixel area (10 micron x 10 micron)
    pd_area = 1e-4                 # photodiode surface area in m^2 from
                                   # Hamamatsu S2281 data sheet
    qe = {}
    wlarrs = {}
    for amp in imutils.allAmps:
        qe[amp] = []
        wlarrs[amp] = []
        for i, wl_nm in enumerate(qe_data.wl):
            try:
                power = qe_data.pd[i]/pd_sph(wl_nm)*ccd_frac(wl_nm)/pd_area
                wl = wl_nm*1e-9
                numerator = 100.*planck*clight*qe_data.medians[amp][i]
                denominator = gains[amp]*qe_data.exptime[i]*power*pixel_area*wl
                if denominator == 0:
                    raise ZeroDivisionError("Zero-valued denominator in "
                                            "qe_value calculation")
                qe_value = numerator/denominator
                qe[amp].append(qe_value)
                wlarrs[amp].append(wl_nm)
            except RuntimeError:
                continue
    for i, amp in enumerate(imutils.allAmps):
        wlarrs[amp] = np.array(wlarrs[amp])
        qe[amp] = np.abs(np.array(qe[amp]))
        plot.curve(wlarrs[amp], qe[amp], oplot=i,
                   xname='wavelength (nm)', yname='QE (%)')
