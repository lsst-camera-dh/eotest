"""
@brief Fe55 K-alpha e-h yield in silicon as a function of temperature.
This function is based on lab measurements from B. G. Lowe &
R. A. Sareen, 2007, NIMA, 576, 367 (L&S).
"""

def fe55_yield(ccdtemp):
    T = 273.2 + ccdtemp  # Convert to Kelvin.
    Ex = 5.889e3         # Mn K-alpha energy (eV)
    slope = -0.000131    # From L&S
    Tref = 270           # L&S figure 4.
    Eref = 3.68
    Epair = slope*(T - Tref) + Eref  # e-h pair creation energy (eV)
    Ne = int(Ex/Epair)
    return Ne, Epair

if __name__ == '__main__':
    import numpy as np
    for ccdtemp in np.linspace(-122.8, 173, 20):
        Ne, Epair = fe55_yield(ccdtemp)
        print "%4i   %4i   %.3f   %4i" % (273.2+ccdtemp, ccdtemp, Epair, Ne)
