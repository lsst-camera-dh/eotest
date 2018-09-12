"""
@brief Fe55 K-alpha e-h yield in silicon as a function of temperature.
This function is based on lab measurements from B. G. Lowe &
R. A. Sareen, 2007, NIMA, 576, 367 (L&S).
"""
from __future__ import print_function


from builtins import object
def pair_energy(ccdtemp):
    """Electron-hole pair creation energy in silicon."""
    T = 273.2 + ccdtemp  # Convert to Kelvin.
    slope = -0.000131    # From L&S given in %/K
    Tref = 270           # L&S figure 4.
    Eref = 3.654
    Eref_sigma = 0.040   # See LSSTTD-141@JIRA
    Epair = (1 + slope*(T - Tref))*Eref  # e-h pair creation energy (eV)
    return Epair, (1 + slope*(T - Tref))*Eref_sigma


class Fe55Yield(object):
    def __init__(self, ccdtemp):
        self.ccdtemp = ccdtemp
        self.pair_energy, self.pair_energy_error = pair_energy(ccdtemp)

    def __call__(self, Ex):
        Ne = Ex/self.pair_energy
        sigma = Ne*self.pair_energy_error/self.pair_energy
        return int(Ne), sigma

    def alpha(self, Ex=5.889e3):
        return self(Ex)

    def beta(self, Ex=6.490e3):
        return self(Ex)


if __name__ == '__main__':
    import numpy as np
    for ccdtemp in np.linspace(-122.8, 173, 20):
        fe55_yield = Fe55Yield(ccdtemp)
        Epair = fe55_yield.pair_energy
        Ne, sigma = fe55_yield.alpha()
        print("%4i   %4i   %.3f   %4i" % (273.2+ccdtemp, ccdtemp, Epair, Ne))
