"""
@brief Generate a standard set of simulated data for a subset of
sensors using a specified directory structure.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os

import numpy as np
import image_utils as imutils
from qe.PhotodiodeResponse import PhotodiodeResponse, CcdIllumination

from sim_tools import *

def mkdir(path):
    try:
        os.makedirs(path)
    except OSError:
        pass

def setup(pars, testtype):
    outputdir = os.path.join(pars.rootdir, pars.sensor_id, testtype, 'data')
    mkdir(outputdir)
    sensor_id = pars.sensor_id.replace('-', '_')
    return outputdir, sensor_id

def simulate_frame(exptime, pars, ccdtemp=-100):
    sensor = CCD(exptime=exptime, gain=pars.system_gain, ccdtemp=ccdtemp,
                 full_well=pars.full_well)
    #
    # Add test-independent effects.
    #
    sensor.add_bias(level=pars.bias_level, sigma=pars.bias_sigma)
    sensor.add_bias(level=0, sigma=pars.read_noise)
    sensor.add_dark_current(pars.dark_current)
    return sensor

def generate_Fe55(pars):
    outputdir, sensor_id = setup(pars, 'xray')
    fe55 = pars.fe55
    for frame in range(fe55.nframes):
        print "Generating Fe55 frame", frame
        #
        # Fe55 exposure
        #
        sensor = simulate_frame(fe55.exptime, pars, ccdtemp=fe55.ccdtemp)
        sensor.add_Fe55_hits(nxrays=fe55.nxrays)
        filename = "%s_fe55_%04is_%03i.fits" % (sensor_id, fe55.exptime, frame)
        sensor.writeto(os.path.join(outputdir, filename))
        #
        # Bias frame
        #
        exptime = 0
        sensor = simulate_frame(exptime, pars, ccdtemp=fe55.ccdtemp)
        filename = "%s_fe55_bias_%03i.fits" % (sensor_id, frame)
        sensor.writeto(os.path.join(outputdir, filename))

def generate_darks(pars):
    outputdir, sensor_id = setup(pars, 'dark')
    darks = pars.darks
    bright_cols = None
    bright_pix = None
    for ccdtemp in darks.ccdtemps:
        #
        # Bias frame
        #
        exptime = 0
        sensor = simulate_frame(exptime, pars)
        filename = "%s_bias_%03i.fits" % (sensor_id, -ccdtemp)
        sensor.writeto(os.path.join(outputdir, filename))
        #
        # Dark frames
        #
        for frame in range(darks.nframes):
            sensor = simulate_frame(darks.exptime, pars)
            if bright_cols is None:
                bright_cols = sensor.generate_bright_cols(darks.bright_ncols)
            if bright_pix is None:
                bright_pix = sensor.generate_bright_pix(darks.bright_npix)
            sensor.add_bright_cols(bright_cols, nsig=darks.bright_nsig)
            sensor.add_bright_pix(bright_pix, nsig=darks.bright_nsig)

            filename = "%s_dark_%03i_%03i.fits" % (sensor_id, -ccdtemp, frame)
            sensor.writeto(os.path.join(outputdir, filename))

def generate_flats(pars):
    outputdir, sensor_id = setup(pars, 'flat')
    flats = pars.flat_fields
    Nes = np.logspace(np.log10(flats.min_charge), np.log10(flats.max_charge),
                      flats.nframes)
    exptimes = np.logspace(np.log10(flats.exptime_min),
                           np.log10(flats.exptime_max),
                           flats.nframes)
    for exptime, Ne in zip(exptimes, Nes):
        intensity = Ne/exptime
        for flat_id in ('flat1', 'flat2'):
            sensor = CCD(exptime=exptime, gain=pars.system_gain,
                         ccdtemp=pars.flat_fields.ccdtemp,
                         full_well=pars.full_well)
            sensor.md['MONDIODE'] = intensity
            sensor.expose_flat(intensity)
            #
            # Test-independent effects. These need to be added after
            # exposure so that full_well can be applied to the
            # collected charge.  (Need to determine if dark current
            # should be applied here or earlier.)
            #
            sensor.add_bias(level=pars.bias_level, sigma=pars.bias_sigma)
            sensor.add_bias(level=0, sigma=pars.read_noise)
            sensor.add_dark_current(pars.dark_current)
            filename = "%(sensor_id)s_flat_%(exptime)06.2fs_%(flat_id)s.fits" \
                       % locals()
            sensor.writeto(os.path.join(outputdir, filename))

def generate_qe_dataset(pars):
    outputdir, sensor_id = setup(pars, 'qe')
    wlscan = pars.wavelength_scan
    pd_sph = PhotodiodeResponse(wlscan.sph_cal_file)
    ccd_frac = CcdIllumination(wlscan.wlscan_file, 
                               wlscan.ccd_cal_file,
                               wlscan.sph_cal_file)
    qe = wlscan.qe
    intensity = wlscan.intensity
    for wl_nm in wlscan.wavelengths:
        sensor = CCD(exptime=wlscan.exptime, gain=pars.system_gain, 
                     ccdtemp=wlscan.ccdtemp)
        sensor.md['MONOWL'] = wl_nm
        #
        # The photodiode current is measured at integrating sphere.
        #
        sensor.md['MONDIODE'] = intensity*pd_sph(wl_nm)/ccd_frac(wl_nm)
        #
        sensor.expose_flat(intensity*qe(wl_nm))
        #
        sensor.add_bias(level=pars.bias_level, sigma=pars.bias_sigma)
        sensor.add_bias(level=0, sigma=pars.read_noise)
        sensor.add_dark_current(pars.dark_current)
        #
        filename = "%(sensor_id)s_qe_%(wl_nm)04.2f.fits" % locals()
        sensor.writeto(os.path.join(outputdir, filename))

if __name__ == '__main__':
    import simulation.sim_params as pars
    try:
        pars.rootdir = os.environ['ROOTDIR']
    except KeyError:
        #pars.rootdir = '/nfs/farm/g/lsst/u1/testData/SIMData'
        pars.rootdir = '/nfs/slac/g/ki/ki18/jchiang/LSST/SensorTests/test_scripts/work'
    pars.sensor_id = '000-00'
#    generate_flats(pars)
#    generate_darks(pars)
#    generate_Fe55(pars)
    generate_qe_dataset(pars)
