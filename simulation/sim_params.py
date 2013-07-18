"""
@brief Example parameter file for simulations.  We may want a more
structured configuration file in the future.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import numpy as np
from qe.PhotodiodeResponse import Interpolator
from sim_tools import xtalk_pattern

class Params(object):
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

bitpix = 16

sensor_id = '000-00'
rootdir = '.'

system_gain = 5.
bias_level = 1e4
bias_sigma = 4
read_noise = 5.
dark_current = 2e-3
full_well = 150000
debug = True

flat_fields = Params(test_type='flat',
                     min_charge=100,
                     max_charge=2e5,
                     exptime_min=1,
                     exptime_max=100,
                     nframes=100,
                     ccdtemp=-95)

traps = Params(test_type='trap',
               Ne=20000,
               exptime=1,
               ccdtemp=-95,
               cycles=100,
               size=200,
               ndefects=100,
               bias_frames=True)

darks = Params(test_type='dark',
               ccdtemp=-95,
               exptime=500,
               nframes=5,
               bright_ncols=1,
               bright_npix=100,
               bright_nsig=5)

fe55 = Params(test_type='fe55',
              nframes=25,
              nxrays=1000,
              exptime=10,
              ccdtemp=-95,
              sigma=0.36)

datapath = lambda x : os.path.join(os.environ['SCRIPTDIR'], 'qe', x)

#qe = lambda wl_nm : 1
qe_curve = np.recfromtxt(datapath('qe_curve.txt')).transpose()
qe = Interpolator(*qe_curve)
wavelength_scan = Params(test_type='lambda',
                         wavelengths=range(400, 1000, 10),
                         exptime=1,
                         ccdtemp=-95,
                         wlscan_file=datapath('WLscan.txt'),
                         ccd_cal_file=datapath('OD142.csv'),
                         sph_cal_file=datapath('OD143.csv'),
                         qe=qe,
                         incident_power=1e-16)
#wavelength_scan.wavelengths.extend((325, 355, 385))
#wavelength_scan.wavelengths.sort()

superflat = Params(test_type='superflat_500',
                   nframes=25,
                   wavelength=500.,
                   exptime=100,
                   pcti=1e-3,
                   scti=1e-3,
                   verbose=True)

spot = Params(test_type='spot',
              exptime=1,
              ccdtemp=-95,
              xtalk_pattern=xtalk_pattern,
              frac_scale=0.02,
              dn=200,
              x=250, y=250, radius=20)

sysnoise = Params(test_type='noise',
                  nframes=10)

sysxtalk = Params(test_type='xtalk',
                  dn=6e4,
                  column=260)
