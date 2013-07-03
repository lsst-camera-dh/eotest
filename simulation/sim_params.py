"""
@brief Example parameter file for simulations.  We may want a more structured
configuration file in the future.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import numpy as np
from qe.PhotodiodeResponse import Interpolator

class Params(object):
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

system_gain = 5.
bias_level = 1e4
bias_sigma = 4
read_noise = 5.
dark_current = 2e-3
full_well = 150000

flat_fields = Params(min_charge=100,
                     max_charge=3e5,
                     exptime_min=1,
                     exptime_max=100,
                     nframes=100,
                     ccdtemp=-100)

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

datapath = lambda x : os.path.join(os.environ['SCRIPTDIR'], 'qe', x)

#qe = lambda wl_nm : 1
qe_curve = np.recfromtxt(datapath('qe_curve.txt')).transpose()
qe = Interpolator(*qe_curve)
wavelength_scan = Params(wavelengths=range(400, 1000, 10),
                         exptime=1,
                         ccdtemp=-95,
                         wlscan_file=datapath('WLscan.txt'),
                         ccd_cal_file=datapath('OD142.csv'),
                         sph_cal_file=datapath('OD143.csv'),
                         qe=qe,
                         incident_power=1e-16)
#wavelength_scan.wavelengths.extend((325, 355, 385))
wavelength_scan.wavelengths.sort()

superflat = Params(nframes=25,
                   exptime=100)

spot = Params(charge_level=9e4)

sysnoise = Params(nframes=10)
