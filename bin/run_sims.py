#!/usr/bin/env python

"""
@brief Generate a standard set of simulated data for a subset of
sensors using a specified directory structure.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
from __future__ import print_function
import os
import copy
import time
import numpy as np
import lsst.eotest.image_utils as imutils
from lsst.eotest.sensor.PhotodiodeResponse import Interpolator
from lsst.eotest.sensor.QE import planck, clight
from lsst.eotest.sensor.sim_tools import *
from lsst.eotest.sensor.ctesim import ctesim


class CcdFactory(object):
    """
    Factory class for creating CCD objects with a common geometry.
    This uses the Borg pattern so that the geometry may be set by
    configuration outside of the calling functions.
    """
    _shared_state = {}

    def __init__(self, geometry=None):
        self.__dict__ = self._shared_state
        if geometry is not None:
            self.geometry = geometry

    def create(self, *args, **kwds):
        my_kwds = copy.deepcopy(kwds)
        try:
            my_kwds['geometry'] = self.geometry
        except AttributeError:
            pass
        return CCD(*args, **my_kwds)


def ccd(*args, **kwds):
    user_factory = CcdFactory()
    return user_factory.create(*args, **kwds)


class AmpIndexDecorator(object):
    def __init__(self, var):
        self.var = var

    def __getitem__(self, i):
        try:
            return self.var[i]
        except TypeError:
            return self.var


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


def simulate_frame(exptime, pars, ccdtemp=-95, set_full_well=True):
    if set_full_well:
        full_well = pars.full_well
    else:
        full_well = None
    sensor = ccd(exptime=exptime, gain=pars.system_gain, ccdtemp=ccdtemp,
                 full_well=full_well)
    #
    # Add test-independent effects.
    #
    sensor.add_bias(level=pars.bias_level, sigma=pars.bias_sigma)
    sensor.add_bias(level=0, sigma=pars.read_noise)
    sensor.add_dark_current(pars.dark_current)
    return sensor


def generate_flats(pars):
    print("Generating flats dataset...")
    flats = pars.flat_fields
    outputdir, sensor_id = setup(pars, flats.test_type)
    Nes = np.logspace(np.log10(flats.min_charge), np.log10(flats.max_charge),
                      flats.nframes)
    exptimes = np.logspace(np.log10(flats.exptime_min),
                           np.log10(flats.exptime_max),
                           flats.nframes)
    for exptime, Ne in zip(exptimes, Nes):
        print("  exptime %06.2fs" % exptime)
        intensity = Ne/exptime
        for flat_id in ('flat1', 'flat2'):
            sensor = ccd(exptime=exptime, gain=pars.system_gain,
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
            filename = ("%s_%s_flat_%06.2fs_%s_%s.fits"
                        % (sensor_id, flats.test_type, exptime, flat_id,
                           time_stamp(debug=pars.debug)))
            sensor.writeto(os.path.join(outputdir, filename),
                           bitpix=pars.bitpix)


def generate_traps(pars):
    print("Generating trap frames...")
    traps = pars.traps
    outputdir, sensor_id = setup(pars, traps.test_type)
    intensity = traps.Ne/traps.exptime
    frames = []
    #
    # First bias image
    #
    frames.append(ccd(exptime=0, gain=pars.system_gain,
                      ccdtemp=traps.ccdtemp))
    #
    # Pocket pumped image
    #
    ppump = ccd(exptime=traps.exptime, gain=pars.system_gain,
                ccdtemp=traps.ccdtemp)
    ppump.expose_flat(intensity)
    ppump.add_traps(traps.ndefects, traps.cycles, traps.size)
    frames.append(ppump)
    #
    # Second bias image
    #
    frames.append(ccd(exptime=0, gain=pars.system_gain,
                      ccdtemp=traps.ccdtemp))
    #
    # Un-pumped image (i.e., a normal flat)
    #
    flat = ccd(exptime=traps.exptime, gain=pars.system_gain,
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
        frame.md['TESTTYPE'] = 'TRAP'
        frame.md['IMGTYPE'] = imgtype.split('_')[0].upper()
        frame.md['LSST_NUM'] = sensor_id
        if imgtype == 'bias_2':
            seqno = '001'
        else:
            seqno = '000'
        filename = ("%s_%s_%s_%s_%s.fits"
                    % (sensor_id, traps.test_type, imgtype.split('_')[0],
                       seqno, time_stamp(debug=pars.debug)))
        frame.writeto(os.path.join(outputdir, filename),
                      bitpix=pars.bitpix)


def generate_darks(pars):
    print("Generating darks dataset...")
    darks = pars.darks
    outputdir, sensor_id = setup(pars, darks.test_type)
    bright_cols = None
    bright_pix = None
    ccdtemp = darks.ccdtemp
    for frame in range(darks.nframes):
        print("  frame", frame)
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
        # Generate bright column and bright pixel locations only once and
        # apply the same sets of locations to each frame.
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
    print("Generating Fe55 dataset...")
    fe55 = pars.fe55
    outputdir, sensor_id = setup(pars, fe55.test_type)
    for frame in range(fe55.nframes):
        print("  frame", frame)
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


class SpherePhotodiodeCurrent(object):
    def __init__(self, pd_ratio_file, pixel_area=1e-10, pd_area=1e-4):
        data = np.recfromtxt(pd_ratio_file, skip_header=1,
                             names='monowl, sens, qe, ccdfrac, foo, bar')
        self.response = Interpolator(data['monowl'],
                                     data['sens']*data['ccdfrac']
                                     * pd_area/pixel_area)

    def __call__(self, incident_power, wl):
        """
        @return Photodiode current (nA)
        @param incident_power Incident power on CCD (W)
        @param wl Wavelength (nm)
        """
        pd_current = incident_power*self.response(wl)*1e9
        return pd_current


def generate_qe_dataset(pars):
    print("Generating wavelength scan dataset...")
    wlscan = pars.wavelength_scan
    outputdir, sensor_id = setup(pars, wlscan.test_type)
    pd_current = SpherePhotodiodeCurrent(wlscan.pd_ratio_file,
                                         pixel_area=wlscan.pixel_area,
                                         pd_area=wlscan.pd_area)
    qe = wlscan.qe
    incident_power = wlscan.incident_power  # J/s per pixel
    for wl_nm in wlscan.wavelengths:
        print("  wavelength %06.1f nm" % wl_nm)
        sensor = ccd(exptime=wlscan.exptime, gain=pars.system_gain,
                     ccdtemp=wlscan.ccdtemp)
        hnu = planck*clight/wl_nm/1e-9    # photon energy (J)
        sensor.md['MONOWL'] = wl_nm
        #
        # The photodiode current is measured at integrating sphere.
        # Units are nA.
        #
        sensor.md['MONDIODE'] = pd_current(incident_power, wl_nm)
        #
        sensor.expose_flat(incident_power*qe(wl_nm)/hnu)
        #
        sensor.add_bias(level=pars.bias_level, sigma=pars.bias_sigma)
        sensor.add_bias(level=0, sigma=pars.read_noise)
        sensor.add_dark_current(pars.dark_current)
        #
        sensor.md['TESTTYPE'] = 'LAMBDA'
        sensor.md['IMGTYPE'] = 'FLAT'
        sensor.md['LSST_NUM'] = sensor_id
        filename = ("%s_%s_flat_%06.1f_%s.fits"
                    % (sensor_id, wlscan.test_type, wl_nm,
                       time_stamp(debug=pars.debug)))
        sensor.writeto(os.path.join(outputdir, filename),
                       bitpix=pars.bitpix)


def generate_superflat(pars):
    print("Generating superflat dataset...")
    superflat = pars.superflat
    outputdir, sensor_id = setup(pars, superflat.test_type)
    tempfile = os.path.join(outputdir, 'superflat_temp.fits')
    dark_cols = None
    dark_pix = None
    #
    # Set incident flux (ph/s/pixel) so that full well is attained in
    # a single exposure (but also disable full well below).
    #
    intensity = pars.full_well*pars.system_gain/superflat.exptime
    for frame in range(superflat.nframes):
        print("  frame", frame)
        sensor = simulate_frame(superflat.exptime, pars, set_full_well=False)
        sensor.expose_flat(intensity)
        # Generate dark column and dark pixel locations only once and
        # apply the same sets of locations to each frame.
        if dark_cols is None:
            dark_cols = sensor.generate_bright_cols(superflat.dark_ncols)
        if dark_pix is None:
            dark_pix = sensor.generate_bright_pix(superflat.dark_npix)
        sensor.set_dark_cols(dark_cols, superflat.dark_frac)
        sensor.set_dark_pix(dark_pix, superflat.dark_frac)
        sensor.md['MONOWL'] = superflat.wavelength
        sensor.md['TESTTYPE'] = 'SFLAT_500'
        sensor.md['IMGTYPE'] = 'FLAT'
        sensor.md['LSST_NUM'] = sensor_id
        sensor.md['PCTI'] = superflat.pcti
        sensor.md['SCTI'] = superflat.scti
        sensor.writeto(tempfile)
        foo = ctesim(tempfile, pcti=superflat.pcti, scti=superflat.scti,
                     verbose=superflat.verbose)
        filename = ("%s_%s_flat_%02i_%s.fits"
                    % (sensor_id, superflat.test_type, frame,
                       time_stamp(debug=pars.debug)))
        foo.writeto(os.path.join(outputdir, filename), overwrite=True)
    #
    # Generate non-uniform illumination correction file.  For these
    # simulations, we set every pixel to unity, i.e., no correction.
    #
    print("  generating the non-uniform illumination correction file...")
    sensor = ccd(exptime=0)
    for amp in imutils.allAmps():
        sensor.segments[amp].image += 1
    filename = ("%s_illumation_correction_%s.fits"
                % (sensor_id, time_stamp(debug=pars.debug)))
    sensor.writeto(os.path.join(outputdir, filename),
                   bitpix=pars.bitpix)


def dark_frame(exptime, ccdtemp, pars):
    frame = ccd(exptime=exptime, gain=pars.system_gain, ccdtemp=ccdtemp)
    frame.add_bias(level=pars.bias_level, sigma=pars.bias_sigma)
    frame.add_bias(level=0, sigma=pars.read_noise)
    frame.add_dark_current(pars.dark_current)
    return frame


def generate_crosstalk_dataset(pars):
    print("Generating spot dataset...")
    spot = pars.spot
    outputdir, sensor_id = setup(pars, spot.test_type)
    if spot.multiaggressor:
        sensors = AmpIndexDecorator(dark_frame(spot.exptime, spot.ccdtemp,
                                               pars))
    else:
        sensors = dict((amp, dark_frame(spot.exptime, spot.ccdtemp, pars))
                       for amp in imutils.allAmps())
    for aggressor in imutils.allAmps():
        print("  aggressor amp", aggressor)
        xtalk_frac = spot.xtalk_pattern(aggressor, spot.frac_scale)
        sensor = sensors[aggressor]
        xx, yy, radius = (AmpIndexDecorator(item)[aggressor] for item
                          in (spot.x, spot.y, spot.radius))
        for amp in sensor.segments:
            if amp == aggressor:
                sensor.segments[amp].add_spot_image(spot.dn, xx, yy, radius)
            else:
                dn = spot.dn*xtalk_frac[amp]
                sensor.segments[amp].add_spot_image(dn, xx, yy, radius)
        #
        sensor.md['TESTTYPE'] = 'SPOT'
        sensor.md['IMGTYPE'] = 'SPOT'
        sensor.md['LSST_NUM'] = sensor_id
    if spot.multiaggressor:
        filename = ("%s_%s_spot_multi_%s.fits" % (sensor_id, spot.test_type,
                                                  time_stamp(debug=pars.debug)))
        sensor.writeto(os.path.join(outputdir, filename), bitpix=pars.bitpix)
    else:
        for aggressor in imutils.allAmps():
            filename = ("%s_%s_spot_%02i_%s.fits" % (sensor_id, spot.test_type,
                                                     aggressor,
                                                     time_stamp(debug=pars.debug)))
            sensors[aggressor].writeto(os.path.join(outputdir, filename),
                                       bitpix=pars.bitpix)


def generate_system_read_noise(pars):
    print("Generating system read noise...")
    sysnoise = pars.sysnoise
    outputdir = system_dir(pars, sysnoise.test_type)
    for frame in range(sysnoise.nframes):
        print("  frame", frame)
        sensor = ccd(exptime=0, gain=pars.system_gain)
        sensor.add_bias(level=pars.bias_level, sigma=pars.bias_sigma)
        filename = ("%s_%02i_%s.fits"
                    % (sysnoise.test_type, frame,
                       time_stamp(debug=pars.debug)))
        sensor.writeto(os.path.join(outputdir, filename),
                       bitpix=pars.bitpix)


def generate_system_crosstalk_dataset(pars):
    print("Generating system crosstalk dataset...")
    sysxtalk = pars.sysxtalk
    outputdir = system_dir(pars, sysxtalk.test_type)
    for aggressor in imutils.allAmps():
        print("  aggressor amp", aggressor)
        sensor = ccd(exptime=0, gain=pars.system_gain)
        sensor.segments[aggressor].add_sys_xtalk_col(sysxtalk.dn,
                                                     sysxtalk.column)
        sensor.add_bias(level=pars.bias_level, sigma=pars.bias_sigma)
        #
        filename = ("%s_%02i_%s.fits"
                    % (sysxtalk.test_type, aggressor,
                       time_stamp(debug=pars.debug)))
        sensor.writeto(os.path.join(outputdir, filename),
                       bitpix=pars.bitpix)


def generate_persistence_dataset(pars):
    print("Generating persistence dataset...")
    persistence = pars.persistence
    outputdir, sensor_id = setup(pars, persistence.test_type)

    # Baseline bias frames
    for frame in range(persistence.num_bias_frames):
        exptime = 0
        sensor = simulate_frame(exptime, pars, ccdtemp=pars.ccdtemp)
        sensor.md['TESTTYPE'] = 'PERSISTENCE'
        sensor.md['IMGTYPE'] = 'BIAS'
        sensor.md['LSST_NUM'] = sensor_id
        filename = ("%s_%s_bias_%02i_%s.fits" %
                    (sensor_id, persistence.test_type, frame,
                     time_stamp(debug=pars.debug)))
        sensor.writeto(os.path.join(outputdir, filename),
                       bitpix=pars.bitpix)

    # Pre-flat darks to measure light leakage
    for dark_frame, exptime in enumerate(persistence.exptimes_presat_darks):
        sensor = simulate_frame(exptime, pars, ccdtemp=pars.ccdtemp)
        sensor.md['TESTTYPE'] = 'PERSISTENCE'
        sensor.md['IMGTYPE'] = 'DARK'
        sensor.md['LSST_NUM'] = sensor_id
        # Include simulated exposure time.
        time.sleep(exptime)
        filename = ("%s_%s_dark_%02i_%s.fits" %
                    (sensor_id, persistence.test_type, dark_frame,
                     time_stamp(debug=pars.debug)))
        sensor.writeto(os.path.join(outputdir, filename), bitpix=pars.bitpix)

    # Saturated flat.
    sensor = simulate_frame(persistence.flat_exptime, pars,
                            ccdtemp=pars.ccdtemp)
    # Set the imaging array in each segment to the full well value.
    for amp in sensor.segments:
        sensor.segments[amp].imarr[:] = persistence.flat_Ne/pars.system_gain
    sensor.md['TESTTYPE'] = 'PERSISTENCE'
    sensor.md['IMGTYPE'] = 'FLAT'
    sensor.md['LSST_NUM'] = sensor_id
    time.sleep(persistence.flat_exptime)
    filename = ("%s_%s_flat_%02i_%s.fits" %
                (sensor_id, persistence.test_type, 0,
                 time_stamp(debug=pars.debug)))
    sensor.writeto(os.path.join(outputdir, filename), bitpix=pars.bitpix)

    # Reference time for exponential decay of deferred charge.
    t0 = utcnow()

    # Post-saturation darks for characterizing persistence.
    for i, exptime in enumerate(persistence.exptimes_postsat_darks):
        frame = dark_frame + i + 1
        sensor = simulate_frame(exptime, pars, ccdtemp=pars.ccdtemp)
        sensor.md['TESTTYPE'] = 'PERSISTENCE'
        sensor.md['IMGTYPE'] = 'DARK'
        sensor.md['LSST_NUM'] = sensor_id
        obs_time = utcnow()
        # Include simulated exposure time.
        time.sleep(exptime)
        # Add deferred charge to darks.
        decay_factor = np.exp(-((obs_time - t0).sec + exptime)
                              / persistence.decay_time)
        for amp in sensor.segments:
            sensor.segments[amp].imarr += \
                persistence.deferred_charge/pars.system_gain*decay_factor
        filename = ("%s_%s_dark_%02i_%s.fits" %
                    (sensor_id, persistence.test_type, frame,
                     time_stamp(debug=pars.debug)))
        sensor.writeto(os.path.join(outputdir, filename),
                       bitpix=pars.bitpix, obs_time=obs_time)


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
        import lsst.eotest.sensor.sim_params as pars
    #
    # Set up the master CcdFactory object with the desired amplifier
    # geometry.
    #
    geometry = AmplifierGeometry(prescan=pars.prescan,
                                 nx=pars.nx,
                                 ny=pars.ny,
                                 detxsize=pars.detxsize,
                                 detysize=pars.detysize)
    master_factory = CcdFactory(geometry=geometry)
    #
    # Generate the desired datasets
    #
    for dataset in pars.datasets:
        exec('generate_%s(pars)' % dataset)
