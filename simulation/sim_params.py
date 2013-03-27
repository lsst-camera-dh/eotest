"""
@brief Example parameter file for simulations.  We may want a more structured
configuration file in the future.

@author J. Chiang <jchiang@slac.stanford.edu>
"""

class Params(object):
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

system_gain = 5.
bias_level = 1e4
bias_sigma = 4
read_noise = 5.
dark_current = 2e-3

flat_fields = Params(min_charge=100,
                     max_charge=2e5,
                     exptime_min=1,
                     exptime_max=100,
                     nframes=50)

pocket_pumping = Params(charge_levels=(100, 200, 500, 1000),
                        bias_frames=True)

darks = Params(ccdtemps=(-100, -95, -90),
               exptime=500,
               nframes=5,
               bright_ncols=1,
               bright_npix=100,
               bright_nsig=5)

fe55 = Params(nframes=25,
              nxrays=1000,
              exptime=10,
              ccdtemp=-95)

wavelength_scan = Params(wavelengths=range(240, 1110, 10))
wavelength_scan.wavelengths.extend((325, 355, 385))
wavelength_scan.wavelengths.sort()

superflat = Params(nframes=25,
                   exptime=100)

spot = Params(charge_level=9e4)

sysnoise = Params(nframes=10)

