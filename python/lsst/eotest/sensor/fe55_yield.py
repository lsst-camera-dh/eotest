"""
@brief Fe55 K-alpha e-h yield in silicon as a function of temperature.
This function is based on lab measurements from B. G. Lowe &
R. A. Sareen, 2007, NIMA, 576, 367 (L&S).
"""

def pair_energy(ccdtemp):
    """Electron-hole pair creation energy in silicon."""
    T = 273.2 + ccdtemp  # Convert to Kelvin.
    slope = -0.000131    # From L&S given in %/K
    Tref = 270           # L&S figure 4.
    Eref = 3.654
    Epair = (1 + slope*(T - Tref))*Eref  # e-h pair creation energy (eV)
    return Epair

class Fe55Yield(object):
    def __init__(self, ccdtemp):
        self.ccdtemp = ccdtemp
        self.pair_energy = pair_energy(ccdtemp)
    def alpha(self, Ex=5.889e3):
        return int(Ex/self.pair_energy)
    def beta(self, Ex=6.490e3):
        return int(Ex/self.pair_energy)

if __name__ == '__main__':
    import numpy as np
    for ccdtemp in np.linspace(-122.8, 173, 20):
        fe55_yield = Fe55Yield(ccdtemp)
        Epair = fe55_yield.pair_energy
        Ne = fe55_yield.alpha()
        print "%4i   %4i   %.3f   %4i" % (273.2+ccdtemp, ccdtemp, Epair, Ne)
