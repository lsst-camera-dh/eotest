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
from simulation.ctesim import ctesim

def time_stamp(gmtime=False, debug=False):
    if debug:
        return "debug"
    if gmtime:
        now = time.gmtime()
    else:
        now = time.localtime()
    time_stamp = ('%04i%02i%02i%02i%02i%02i'
                  % (now.tm_year, now.tm_mon, now.tm_mday,
                     now.tm_hour, now.tm_min, now.tm_sec))
    return time_stamp

def date_stamp(gmtime=False, debug=False):
    if debug:
        return "debug"
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

def system_dir(pars, testtype):
    outputdir = os.path.join(pars.rootdir, 'system', testtype,
                             date_stamp(debug=pars.debug))
    mkdir(outputdir)
    return outputdir

def setup(pars, testtype):
    outputdir = os.path.join(pars.rootdir, 'sensorData', pars.sensor_id,
                             testtype, date_stamp(debug=pars.debug))
    mkdir(outputdir)
    sensor_id = pars.sensor_id
    return outputdir, sensor_id

def simulate_frame(exptime, pars, ccdtemp=-100, set_full_well=True):
    if set_full_well:
        full_well = pars.full_well
    else:
        full_well = None
    sensor = CCD(exptime=exptime, gain=pars.system_gain, ccdtemp=ccdtemp,
                 full_well=full_well)
    #
    # Add test-independent effects.
    #
    sensor.add_bias(level=pars.bias_level, sigma=pars.bias_sigma)
    sensor.add_bias(level=0, sigma=pars.read_noise)
    sensor.add_dark_current(pars.dark_current)
    return sensor

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
        print "  exptime %06.2fs" % exptime
        intensity = Ne/exptime
        for flat_id in ('flat1', 'flat2'):
            sensor = CCD(exptime=exptime, gain=pars.system_gain,
                         ccdtemp=pars.flat_fields.ccdtemp,
                         full_well=pars.full_well)
            sensor.md['MONDIODE'] = intensity
            sensor.md['TESTTYPE'] = 'FLAT'
            sensor.md['IMGTYPE'] = 'FLAT'
            sensor.md['LSST_NUM'] = sensor_id
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
            filename = ("%s_%s_%06.2fs_%s_%s.fits" 
                        % (sensor_id, flats.test_type, exptime, flat_id,
                           time_stamp(debug=pars.debug)))
            sensor.writeto(os.path.join(outputdir, filename),
                           bitpix=pars.bitpix)

def generate_traps(pars):
    print "Generating trap frames..."
    traps = pars.traps
    outputdir, sensor_id = setup(pars, traps.test_type)
    intensity = traps.Ne/traps.exptime
    frames = []
    #
    # First bias image
    #
    frames.append(CCD(exptime=0, gain=pars.system_gain,
                      ccdtemp=traps.ccdtemp))
    #
    # Pocket pumped image
    #
    ppump = CCD(exptime=traps.exptime, gain=pars.system_gain,
                ccdtemp=traps.ccdtemp)
    ppump.expose_flat(intensity)
    ppump.add_traps(traps.ndefects, traps.cycles, traps.size)
    frames.append(ppump)
    #
    # Second bias image
    #
    frames.append(CCD(exptime=0, gain=pars.system_gain,
                      ccdtemp=traps.ccdtemp))
    #
    # Un-pumped image (i.e., a normal flat)
    #
    flat = CCD(exptime=traps.exptime, gain=pars.system_gain,
               ccdtemp=traps.ccdtemp)
    flat.expose_flat(intensity)
    frames.append(flat)
    #
    # Add bias and dark current (accounting for zero exposure in bias
    # frames) and write.
    #
    for frame, imgtype in zip(frames, ('bias_1', 'ppump', 'bias_2', 'flat')):
        frame.add_bias(level=pars.bias_level, sigma=pars.bias_sigma)
        frame.add_bias(level=0, sigma=pars.read_noise)
        frame.add_dark_current(pars.dark_current)
        frame.md['TESTTYPE'] = 'PPUMP'
        frame.md['IMGTYPE'] = imgtype.split('_')[0].upper()
        frame.md['LSST_NUM'] = sensor_id
        filename = ("%s_%s_%s_%s.fits" 
                    % (sensor_id, traps.test_type, imgtype,
                       time_stamp(debug=pars.debug)))
        frame.writeto(os.path.join(outputdir, filename),
                      bitpix=pars.bitpix)

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
        sensor.md['TESTTYPE'] = 'DARK'
        sensor.md['IMGTYPE'] = 'BIAS'
        sensor.md['LSST_NUM'] = sensor_id
        filename = ("%s_%s_bias_%02i_%s.fits"
                    % (sensor_id, darks.test_type, frame,
                       time_stamp(debug=pars.debug)))
        sensor.writeto(os.path.join(outputdir, filename),
                       bitpix=pars.bitpix)
        #
        # Dark frame
        #
        sensor = simulate_frame(darks.exptime, pars)
        #
        # Compute desired signal level in terms of image stdev.
        #
        sigma = sensor.segments[1].sigma()
        nsig = darks.bright_Ne_per_sec*darks.exptime/pars.system_gain/sigma
        if bright_cols is None:
            bright_cols = sensor.generate_bright_cols(darks.bright_ncols)
        if bright_pix is None:
            bright_pix = sensor.generate_bright_pix(darks.bright_npix)
        sensor.add_bright_cols(bright_cols, nsig=nsig)
        sensor.add_bright_pix(bright_pix, nsig=nsig)
        sensor.md['TESTTYPE'] = 'DARK'
        sensor.md['IMGTYPE'] = 'DARK'
        sensor.md['LSST_NUM'] = sensor_id
        filename = ("%s_%s_dark_%02i_%s.fits"
                    % (sensor_id, darks.test_type, frame,
                       time_stamp(debug=pars.debug)))
        sensor.writeto(os.path.join(outputdir, filename),
                       bitpix=pars.bitpix)

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
        sensor.add_Fe55_hits(nxrays=fe55.nxrays, sigma=fe55.sigma)
        sensor.md['TESTTYPE'] = 'FE55'
        sensor.md['IMGTYPE'] = 'FE55'
        sensor.md['LSST_NUM'] = sensor_id
        filename = ("%s_%s_fe55_%02i_%s.fits"
                    % (sensor_id, fe55.test_type, frame,
                       time_stamp(debug=pars.debug)))
        sensor.writeto(os.path.join(outputdir, filename),
                           bitpix=pars.bitpix)
        #
        # Bias frame
        #
        exptime = 0
        sensor = simulate_frame(exptime, pars, ccdtemp=fe55.ccdtemp)
        sensor.md['TESTTYPE'] = 'FE55'
        sensor.md['IMGTYPE'] = 'BIAS'
        sensor.md['LSST_NUM'] = sensor_id
        filename = ("%s_%s_bias_%02i_%s.fits"
                    % (sensor_id, fe55.test_type, frame,
                       time_stamp(debug=pars.debug)))
        sensor.writeto(os.path.join(outputdir, filename),
                           bitpix=pars.bitpix)

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
        print "  wavelength %06.1f nm" % wl_nm
        sensor = CCD(exptime=wlscan.exptime, gain=pars.system_gain, 
                     ccdtemp=wlscan.ccdtemp)
        hnu = planck*clight/wl_nm/1e-9    # photon energy (J)
        sensor.md['MONOWL'] = wl_nm
        #
        # The photodiode current is measured at integrating sphere.
        # Units are nA.
        #
        sensor.md['MONDIODE'] = (incident_power*pd_sph(wl_nm)/ccd_frac(wl_nm)
                                 *(pars.pd_area/pars.pixel_area))*1e9
        #
        sensor.expose_flat(incident_power*qe(wl_nm)/hnu)
        #
        sensor.add_bias(level=pars.bias_level, sigma=pars.bias_sigma)
        sensor.add_bias(level=0, sigma=pars.read_noise)
        sensor.add_dark_current(pars.dark_current)
        #
        sensor.md['TESTTYPE'] = 'QE'
        sensor.md['IMGTYPE'] = 'FLAT'
        sensor.md['LSST_NUM'] = sensor_id
        filename = ("%s_%s_%06.1f_%s.fits"
                    % (sensor_id, wlscan.test_type, wl_nm,
                       time_stamp(debug=pars.debug)))
        sensor.writeto(os.path.join(outputdir, filename),
                       bitpix=pars.bitpix)

def generate_superflat(pars):
    print "Generating superflat dataset..."
    superflat = pars.superflat
    outputdir, sensor_id = setup(pars, superflat.test_type)
    tempfile = os.path.join(outputdir, 'superflat_temp.fits')
    #
    # Set incident flux (ph/s/pixel) so that full well is attained in
    # a single exposure (but also disable full well below).
    #
    intensity = pars.full_well*pars.system_gain/superflat.exptime
    for frame in range(superflat.nframes):
        print "  frame", frame
        sensor = simulate_frame(superflat.exptime, pars, set_full_well=False)
        sensor.expose_flat(intensity)
        sensor.md['MONOWL'] = superflat.wavelength
        sensor.md['TESTTYPE'] = 'SFLAT'
        sensor.md['IMGTYPE'] = 'FLAT'
        sensor.md['LSST_NUM'] = sensor_id
        sensor.writeto(tempfile)
        foo = ctesim(tempfile, pcti=superflat.pcti, scti=superflat.scti,
                     verbose=superflat.verbose)
        filename = ("%s_%s_%02i_%s.fits" 
                    % (sensor_id, superflat.test_type, frame,
                       time_stamp(debug=pars.debug)))
        foo.writeto(os.path.join(outputdir, filename), clobber=True)
    #
    # Generate non-uniform illumination correction file.  For these
    # simulations, we set every pixel to unity, i.e., no correction.
    #
    print "  generating the non-uniform illumination correction file..."
    sensor = CCD(exptime=0)
    for amp in imutils.allAmps:
        sensor.segments[amp].image += 1
    filename = ("%s_illumation_correction_%s.fits"
                % (sensor_id, time_stamp(debug=pars.debug)))
    sensor.writeto(os.path.join(outputdir, filename),
                   bitpix=pars.bitpix)

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
        sensor.md['TESTTYPE'] = 'SPOT'
        sensor.md['IMGTYPE'] = 'SPOT'
        sensor.md['LSST_NUM'] = sensor_id
        filename = ("%s_%s_%02i_%s.fits"
                    % (sensor_id, spot.test_type, aggressor,
                       time_stamp(debug=pars.debug)))
        sensor.writeto(os.path.join(outputdir, filename),
                       bitpix=pars.bitpix)

def generate_system_read_noise(pars):
    print "Generating system read noise..."
    sysnoise = pars.sysnoise
    outputdir = system_dir(pars, sysnoise.test_type)
    for frame in range(sysnoise.nframes):
        print "  frame", frame
        sensor = CCD(exptime=0, gain=pars.system_gain)
        sensor.add_bias(level=pars.bias_level, sigma=pars.bias_sigma)
        filename = ("%s_%02i_%s.fits"
                    % (sysnoise.test_type, frame,
                       time_stamp(debug=pars.debug)))
        sensor.writeto(os.path.join(outputdir, filename),
                       bitpix=pars.bitpix)

def generate_system_crosstalk_dataset(pars):
    print "Generating system crosstalk dataset..."
    sysxtalk = pars.sysxtalk
    outputdir = system_dir(pars, sysxtalk.test_type)
    for aggressor in imutils.allAmps:
        print "  aggressor amp", aggressor
        sensor = CCD(exptime=0, gain=pars.system_gain)
        sensor.segments[aggressor].add_sys_xtalk_col(sysxtalk.dn,
                                                     sysxtalk.column)
        sensor.add_bias(level=pars.bias_level, sigma=pars.bias_sigma)
        #
        filename = ("%s_%02i_%s.fits"
                    % (sysxtalk.test_type, aggressor,
                       time_stamp(debug=pars.debug)))
        sensor.writeto(os.path.join(outputdir, filename),
                       bitpix=pars.bitpix)

if __name__ == '__main__':
    import sys
    #
    # Import the parameter file.
    #
    if sys.argv[1:]:
        parfile = sys.argv[1]
        tmp_path = os.path.join(os.path.split(parfile)[:-1])[0]
        sys.path.insert(0, tmp_path)
        exec('import %s as pars' % os.path.basename(parfile).strip('.py'))
    else:
        import simulation.sim_params as pars
    #
    # Generate the desired datasets
    #
    for dataset in pars.datasets:
        exec('generate_%s(pars)' % dataset)
