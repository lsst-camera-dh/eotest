"""
@brief Classes to encapsulate NIST calibrations for photodiodes and
wavelength scans of photodiode readings at integrating sphere and CCD
location.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import bisect
import numpy as np

class Interpolator(object):
    def __init__(self, xx, yy):
        self.xx = xx
        self.yy = yy
        dx = xx[1:] - xx[:-1]
        if max(dx) == min(dx):
            self.dx = dx[0]
        else:
            self.dx = None
    def __call__(self, x):
        try:
            return np.array([self._value(xval) for xval in x])
        except TypeError:
            return self._value(x)
    def _value(self, x):
        if x < self.xx[0] or x > self.xx[-1]:
            raise RuntimeError("Requested x-value is out-of-range.")
        elif x == self.xx[0]:
            return self.yy[0]
        elif x == self.xx[-1]:
            return self.yy[-1]
        
        if self.dx is not None:
            indx = int((x - self.xx[0])/self.dx)
        else:
            indx = bisect.bisect(self.xx, x) - 1
        indx = max(0, min(len(self.xx)-2, indx))
        value = ( (x - self.xx[indx])/
                  (self.xx[indx+1] - self.xx[indx])*
                  (self.yy[indx+1] - self.yy[indx]) + self.yy[indx] )
        return value

class PhotodiodeResponse(Interpolator):
    def __init__(self, calfile, skip_header=5,):
        data = np.recfromtxt(calfile, delimiter=',', skip_header=skip_header,
                             names="wavelength, resp, pct_error, none, none2")
        self.wavelength = data['wavelength']
        self.resp = data['resp']
        Interpolator.__init__(self, self.wavelength, self.resp)

class CcdIllumination(Interpolator):
    """
    This functor returns the ratio of the intensity at the photodiode
    at the CCD position to the intensity at the photodiode located on
    the integrating sphere as a function of wavelength.
    """
    def __init__(self, scan_file, ccd_cal_file, sph_cal_file):
        data = np.recfromtxt(scan_file, names=True)
        self.wavelengths = data['wl']
        pd_ccd = PhotodiodeResponse(ccd_cal_file)
        pd_sph = PhotodiodeResponse(sph_cal_file)
        intsphere = data['intsphere']/pd_sph(self.wavelengths)
        ccdpos = data['ccdpos']/pd_ccd(self.wavelengths)
        self.ccdfrac = ccdpos/intsphere
        Interpolator.__init__(self, self.wavelengths, self.ccdfrac)

if __name__ == '__main__':
    import pylab_plotter as plot

    pd1 = PhotodiodeResponse('OD142.csv')
    pd2 = PhotodiodeResponse('OD143.csv')

    lamb = np.linspace(400, 1000, 200)
    plot.curve(lamb, pd1(lamb), xname='wavelength (nm)',
               yname='Photodiode response')
    plot.curve(lamb, [pd2(x) for x in lamb], oplot=1, color='r')

    ccd_frac = CcdIllumination('WLscan.txt', 'OD142.csv', 'OD143.csv')

    plot.curve(lamb, ccd_frac(lamb), xname='wavelength (nm)',
               yname='ratio of illumination at CCD to integrating sphere')
