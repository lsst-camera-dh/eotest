"""
@brief Example parameter file for simulations.  We may want a more
structured configuration file in the future.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
from __future__ import absolute_import
import os
import numpy as np
from .PhotodiodeResponse import Interpolator
from .multiaggressor_tools import multiaggressor_amplifier_coords
from .sim_tools import CrosstalkPattern


class Params(object):
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def __setattr__(self, key, value):
        """
        Perform validation by checking that a parameter name already
        exists before trying to set it.
        """
        if not hasattr(self, key):
            raise AttributeError("Cannot add a new attribute to 'Params' object.")
        self.__dict__[key] = value


pd_area = 1e-4          # Sensitive area of photodiode
pixel_area = 1e-10      # Nominal pixel area (10 micron x 10 micron)

bitpix = 16             # Only bitpix = 16, -32 are supported.

sensor_id = '000-00'
rootdir = '.'

system_gain = 5.
bias_level = 1e4
bias_sigma = 4
read_noise = 5.
dark_current = 2e-3
full_well = 150000
ccdtemp = -95

#
# Geometry of amplifiers in a sensor.  These are nomimal values for
# e2v devices.
#
prescan = 10
nx = 512
ny = 2002
detxsize = 4336
detysize = 4044

#
# If debug=true, then 'debug' will be used in place of the date and time
# stamp strings in the directory and file names.
#
debug = True

#
# Leave desired datasets uncommented to generate them.
#
datasets = [
    'flats',
    'traps',
    'darks',
    'Fe55',
    'qe_dataset',
    'superflat',
    'crosstalk_dataset',
    'system_read_noise',
    'system_crosstalk_dataset',
    'persistence_dataset'
]

flat_fields = Params(test_type='flat',
                     min_charge=100,
                     max_charge=2e5,
                     exptime_min=1,
                     exptime_max=100,
                     nframes=100,
                     ccdtemp=ccdtemp)

traps = Params(test_type='trap',
               Ne=20000,
               exptime=1,
               ccdtemp=ccdtemp,
               cycles=100,
               size=200,
               ndefects=100,
               bias_frames=True)

darks = Params(test_type='dark',
               ccdtemp=ccdtemp,
               exptime=500,
               nframes=5,
               bright_ncols=1,
               bright_npix=100,
               bright_Ne_per_sec=10)

fe55 = Params(test_type='fe55',
              nframes=25,
              nxrays=1000,
              exptime=10,
              ccdtemp=ccdtemp,
              sigma=0.36)
#
# Determine the path to the qe subdirectory
#
_qe_dir = os.path.join(os.environ['EOTEST_DIR'], 'data', 'qe')


def qe_path(x): return os.path.join(_qe_dir, x)


#qe = lambda wl_nm : 1   # 100% quantum efficiency
qe_curve = np.recfromtxt(qe_path('sim/qe_curve.txt')).transpose()
qe = Interpolator(*qe_curve)
wavelength_scan = Params(test_type='lambda',
                         wavelengths=range(320, 1110, 10),
                         exptime=1,
                         ccdtemp=ccdtemp,
                         pd_ratio_file=qe_path('BNL/pd_Cal_mar2013_v1.txt'),
                         qe=qe,
                         incident_power=1e-16,
                         pixel_area=1e-10,
                         pd_area=1e-4)
wavelength_scan.wavelengths.extend((325, 355, 385))
wavelength_scan.wavelengths.sort()

superflat = Params(test_type='sflat_500',
                   nframes=25,
                   wavelength=500.,
                   exptime=100,
                   pcti=1e-3,
                   scti=1e-3,
                   dark_npix=100,
                   dark_ncols=1,
                   dark_frac=0.5,
                   verbose=True)

xpos, ypos = multiaggressor_amplifier_coords(nx, ny)

spot = Params(test_type='spot',
              exptime=1,
              ccdtemp=ccdtemp,
              xtalk_pattern=CrosstalkPattern(),
              frac_scale=0.02,
              dn=200,
              radius=20,
              #              x=250, y=250,
              #              multiaggressor=False,
              x=xpos, y=ypos,
              multiaggressor=True,
              )

sysnoise = Params(test_type='noise',
                  nframes=10)

sysxtalk = Params(test_type='xtalk',
                  dn=6e4,
                  column=260)

persistence = Params(test_type='persistence',
                     num_bias_frames=5,
                     exptimes_presat_darks=(20, 20, 20),
                     flat_exptime=40,
                     flat_Ne=full_well,
                     exptimes_postsat_darks=(0, 0, 0, 5, 5, 5, 10, 10, 20, 20),
                     deferred_charge=100,  # e-/pixel
                     decay_time=1e2)       # seconds
