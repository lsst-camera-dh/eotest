"""
@brief Generate a standard set of simulated data for a subset of
sensors using a specified directory structure.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import time
import numpy as np
import image_utils as imutils
from qe.PhotodiodeResponse import PhotodiodeResponse, CcdIllumination
from qe.QE import planck, clight
from simulation.sim_tools import *

pd_area = 1e-4          # Sensitive area of photodiode
pixel_area = 1e-10      # Nominal pixel area (10 micron x 10 micron)

def time_stamp(gmtime=False):
    if gmtime:
        now = time.gmtime()
    else:
        now = time.localtime()
    time_stamp = ('%04i%02i%02i%02i%02i%02i'
                  % (now.tm_year, now.tm_mon, now.tm_mday,
                     now.tm_hour, now.tm_min, now.tm_sec))
    return time_stamp

def date_stamp(gmtime=False):
    if gmtime:
        now = time.gmtime()
    else:
        now = time.localtime()
    date_stamp = ('%02i%02i%02i-%02i%02i%02i'
                  % (now.tm_year % 100, now.tm_mon, now.tm_mday,
                     now.tm_hour, now.tm_min, now.tm_sec))
    return date_stamp

def mkdir(path):
    try:
        os.makedirs(path)
    except OSError:
        pass

def setup(pars, testtype):
    outputdir = os.path.join(pars.rootdir, 'sensorData', pars.sensor_id,
                             testtype, date_stamp())
    mkdir(outputdir)
    sensor_id = pars.sensor_id
    return outputdir, sensor_id

def simulate_frame(exptime, pars, ccdtemp=-100, set_full_well=True):
    if set_full_well:
        full_well = pars.full_well
    sensor = CCD(exptime=exptime, gain=pars.system_gain, ccdtemp=ccdtemp,
                 full_well=full_well)
    #
    # Add test-independent effects.
    #
    sensor.add_bias(level=pars.bias_level, sigma=pars.bias_sigma)
    sensor.add_bias(level=0, sigma=pars.read_noise)
    sensor.add_dark_current(pars.dark_current)
    return sensor

def generate_Fe55(pars):
    print "Generating Fe55 dataset..."
    fe55 = pars.fe55
    outputdir, sensor_id = setup(pars, fe55.test_type)
    for frame in range(fe55.nframes):
        print "  frame", frame
        #
        # Fe55 exposure
        #
        sensor = simulate_frame(fe55.exptime, pars, ccdtemp=fe55.ccdtemp)
        sensor.add_Fe55_hits(nxrays=fe55.nxrays)
        filename = ("%s_%s_fe55_%02i_%s.fits"
                    % (sensor_id, fe55.test_type, frame, time_stamp()))
        sensor.writeto(os.path.join(outputdir, filename))
        #
        # Bias frame
        #
        exptime = 0
        sensor = simulate_frame(exptime, pars, ccdtemp=fe55.ccdtemp)
        filename = ("%s_%s_bias_%02i_%s.fits"
                    % (sensor_id, fe55.test_type, frame, time_stamp()))
        sensor.writeto(os.path.join(outputdir, filename))

def generate_darks(pars):
    print "Generating darks dataset..."
    darks = pars.darks
    outputdir, sensor_id = setup(pars, darks.test_type)
    bright_cols = None
    bright_pix = None
    ccdtemp = darks.ccdtemp
    for frame in range(darks.nframes):
        print "  frame", frame
        #
        # Bias frame
        #
        exptime = 0
        sensor = simulate_frame(exptime, pars)
        filename = ("%s_%s_bias_%02i_%s.fits"
                    % (sensor_id, darks.test_type, frame, time_stamp()))
        sensor.writeto(os.path.join(outputdir, filename))
        #
        # Dark frame
        #
        sensor = simulate_frame(darks.exptime, pars)
        if bright_cols is None:
            bright_cols = sensor.generate_bright_cols(darks.bright_ncols)
        if bright_pix is None:
            bright_pix = sensor.generate_bright_pix(darks.bright_npix)
        sensor.add_bright_cols(bright_cols, nsig=darks.bright_nsig)
        sensor.add_bright_pix(bright_pix, nsig=darks.bright_nsig)
        filename = ("%s_%s_dark_%02i_%s.fits"
                    % (sensor_id, darks.test_type, frame, time_stamp()))
        sensor.writeto(os.path.join(outputdir, filename))

def generate_flats(pars):
    print "Generating flats dataset..."
    flats = pars.flat_fields
    outputdir, sensor_id = setup(pars, flats.test_type)
    Nes = np.logspace(np.log10(flats.min_charge), np.log10(flats.max_charge),
                      flats.nframes)
    exptimes = np.logspace(np.log10(flats.exptime_min),
                           np.log10(flats.exptime_max),
                           flats.nframes)
    for exptime, Ne in zip(exptimes, Nes):
        print "  exptime %06.2fs", exptime
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
            filename = ("%s_%s_%06.2fs_%s_%s.fits" \
                       % (sensor_id, flats.test_type, exptime, flat_id,
                          time_stamp()))
            sensor.writeto(os.path.join(outputdir, filename))

def generate_qe_dataset(pars):
    print "Generating wavelength scan dataset..."
    wlscan = pars.wavelength_scan
    outputdir, sensor_id = setup(pars, wlscan.test_type)
    pd_sph = PhotodiodeResponse(wlscan.sph_cal_file)
    ccd_frac = CcdIllumination(wlscan.wlscan_file, 
                               wlscan.ccd_cal_file,
                               wlscan.sph_cal_file)
    qe = wlscan.qe
    incident_power = wlscan.incident_power  # J/s per pixel
    for wl_nm in wlscan.wavelengths:
        print "  wavelength %06.1 nm" % wl_nm
        sensor = CCD(exptime=wlscan.exptime, gain=pars.system_gain, 
                     ccdtemp=wlscan.ccdtemp)
        hnu = planck*clight/wl_nm/1e-9    # photon energy (J)
        sensor.md['MONOWL'] = wl_nm
        #
        # The photodiode current is measured at integrating sphere.
        # Units are nA.
        #
        sensor.md['MONDIODE'] = (incident_power*pd_sph(wl_nm)/ccd_frac(wl_nm)
                                 *(pd_area/pixel_area))*1e9
        #
        sensor.expose_flat(incident_power*qe(wl_nm)/hnu)
        #
        sensor.add_bias(level=pars.bias_level, sigma=pars.bias_sigma)
        sensor.add_bias(level=0, sigma=pars.read_noise)
        sensor.add_dark_current(pars.dark_current)
        #
        filename = ("%s_%s_%06.1f_%s.fits"
                    % (sensor_id, wlscan.test_type, wl_nm, time_stamp()))
        sensor.writeto(os.path.join(outputdir, filename))

def generate_crosstalk_dataset(pars):
    print "Generating spot dataset..."
    spot = pars.spot
    outputdir, sensor_id = setup(pars, spot.test_type)
    for aggressor in imutils.allAmps:
        print "  aggressor amp", aggressor
        sensor = CCD(exptime=spot.exptime, gain=pars.system_gain,
                     ccdtemp=spot.ccdtemp)
        xtalk_frac = spot.xtalk_pattern(aggressor, spot.frac_scale)
        for amp in sensor.segments:
            if amp == aggressor:
                sensor.segments[amp].add_spot_image(spot.dn, spot.x, spot.y,
                                                    spot.radius)
            else:
                dn = spot.dn*xtalk_frac[amp]
                sensor.segments[amp].add_spot_image(dn, spot.x, spot.y,
                                                    spot.radius)
        sensor.add_bias(level=pars.bias_level, sigma=pars.bias_sigma)
        sensor.add_bias(level=0, sigma=pars.read_noise)
        sensor.add_dark_current(pars.dark_current)
        #
        filename = ("%s_%s_%02i_%s.fits"
                    % (sensor_id, spot.test_type, aggressor, time_stamp()))
        sensor.writeto(os.path.join(outputdir, filename))

if __name__ == '__main__':
    import simulation.sim_params as pars
    try:
        pars.rootdir = os.environ['ROOTDIR']
    except KeyError:
        #pars.rootdir = '/nfs/farm/g/lsst/u1/testData/SIMData'
        pars.rootdir = '/nfs/slac/g/ki/ki18/jchiang/LSST/SensorTests/test_scripts/work'
    pars.sensor_id = '000-00'
    generate_flats(pars)
    generate_darks(pars)
    generate_Fe55(pars)
    generate_qe_dataset(pars)
    generate_crosstalk_dataset(pars)
