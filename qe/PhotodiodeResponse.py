"""
@brief Class to encapsulate NIST calibrations for photodiodes.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import bisect
import numpy as np

class PhotodiodeResponse(object):
    def __init__(self, calfile, skip_header=5,):
        data = np.recfromtxt(calfile, delimiter=',', skip_header=skip_header,
                             names="wavelength, resp, pct_error, none, none2")
        self.wavelength = data['wavelength']
        self.resp = data['resp']
        self.pct_error = data['pct_error']
        dwl = self.wavelength[1:] - self.wavelength[:-1]
        if max(dwl) == min(dwl):
            self.dwl = dwl[0]
        else:
            self.dwl = None
    def __call__(self, wl):
        if wl < self.wavelength[0] or wl > self.wavelength[-1]:
            raise RuntimeError("Requested wavelength is outside " +
                               "of calibrated range.")
        elif wl == self.wavelength[0]:
            return self.resp[0]
        elif wl == self.wavelength[-1]:
            return self.resp[-1]
        
        if self.dwl is not None:
            indx = int((wl - self.wavelength[0])/self.dwl)
        else:
            indx = bisect.bisect(self.wavelength, wl) - 1
        indx = max(0, min(len(self.wavelength)-2, indx))
        return ( (wl - self.wavelength[indx])/
                 (self.wavelength[indx+1] - self.wavelength[indx])*
                 (self.resp[indx+1] - self.resp[indx]) + self.resp[indx] )

if __name__ == '__main__':
    import pylab_plotter as plot

    pd1 = PhotodiodeResponse('OD142.csv')

    lamb = np.linspace(350, 1100, 200)
    plot.curve(lamb, [pd1(x) for x in lamb], xname='wavelength (nm)',
               yname='Photodiode response')
    
