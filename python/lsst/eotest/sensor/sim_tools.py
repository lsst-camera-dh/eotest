"""
@brief Tools to create simulated CCD segment exposures under ideal
conditions.  Darks, flats, Fe55, etc..
"""
from __future__ import print_function
from __future__ import absolute_import
import os
from collections import OrderedDict
import datetime
import astropy.time

import numpy as np
import numpy.random as random
import astropy.io.fits as fits
from lsst.eotest.fitsTools import fitsWriteto

import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.display.ds9 as ds9

import lsst.eotest.image_utils as imutils
from .fe55_yield import Fe55Yield
from .fits_headers import fits_headers
from .AmplifierGeometry import AmplifierGeometry
from .crosstalk import CrosstalkMatrix

_sqrt2 = np.sqrt(2.)


def utcnow(dt=0):
    now = datetime.datetime.now() + datetime.timedelta(seconds=dt)
    return astropy.time.Time(now.isoformat(), format='isot', scale='utc')


class CrosstalkPattern(object):
    def __init__(self, infile=None):
        if infile is not None:
            print("Using %s for the crosstalk pattern" % infile)
            xtalk_matrix = CrosstalkMatrix(infile)
            self.matrix = xtalk_matrix.matrix
        else:
            # Use default matrix.
            pattern = (0.01, 0.02, 1, 0.02, 0.01)
            offsets = (-2, -1, 0, 1, 2)
            namps = len(imutils.allAmps())
            self.matrix = np.zeros((namps, namps), dtype=np.float)
            for agg in imutils.allAmps():
                for offset, value in zip(offsets, pattern):
                    vic = agg + offset
                    if (vic in range(1, 9) and agg in range(1, 9) or
                            vic in range(9, 17) and agg in range(9, 17)):
                        self.matrix[agg-1][vic-1] = value

    def __call__(self, aggressor, frac_scale=None):
        return dict([(amp, victim) for amp, victim
                     in zip(imutils.allAmps(), self.matrix[aggressor-1])])


def xtalk_pattern(aggressor, frac_scale=0.02):
    xtalk_frac = {}
    nside = len(imutils.allAmps())/2
    for victim in imutils.allAmps():
        if (victim != aggressor and
                (victim-1)/nside == (aggressor-1)/nside):
            dist = abs(victim - aggressor)
            xtalk_frac[victim] = frac_scale/dist**2
        else:
            xtalk_frac[victim] = 0
    return xtalk_frac


class CCD(object):
    dtypes = dict([(-32, np.float32), (16, np.int16)])

    def __init__(self, exptime=1, gain=5, ccdtemp=-95, full_well=None,
                 geometry=AmplifierGeometry(), amps=None):
        self.segments = OrderedDict()
        if amps is None:
            amps = imutils.allAmps()
        for amp in amps:
            self.segments[amp] = SegmentExposure(exptime=exptime,
                                                 gain=gain,
                                                 ccdtemp=ccdtemp,
                                                 full_well=full_well,
                                                 geometry=geometry)
        self.md = dict()

    def add_bias(self, level=1e4, sigma=4):
        for amp in self.segments:
            self.segments[amp].add_bias(level=level, sigma=sigma)

    def add_dark_current(self, level=2e-3):
        for amp in self.segments:
            self.segments[amp].add_dark_current(level=level)

    def expose_flat(self, intensity=0):
        for amp in self.segments:
            self.segments[amp].expose_flat(intensity=intensity)

    def add_Fe55_hits(self, nxrays=1000, beta_frac=0.12, sigma=None):
        for amp in self.segments:
            self.segments[amp].add_Fe55_hits(nxrays=nxrays,
                                             beta_frac=beta_frac,
                                             sigma=sigma)

    def generate_bright_cols(self, ncols=1):
        bright_cols = OrderedDict()
        for amp in self.segments:
            bright_cols[amp] = self.segments[amp].generate_bright_cols(ncols)
        return bright_cols

    def set_dark_cols(self, dark_cols, frac):
        for amp in self.segments:
            self.segments[amp].set_dark_cols(dark_cols[amp], frac)

    def add_bright_cols(self, bright_cols, nsig=5):
        for amp in self.segments:
            self.segments[amp].add_bright_cols(bright_cols[amp], nsig=nsig)

    def generate_bright_pix(self, npix=100):
        bright_pix = OrderedDict()
        for amp in self.segments:
            bright_pix[amp] = self.segments[amp].generate_bright_pix(npix)
        return bright_pix

    def set_dark_pix(self, dark_pix, frac):
        for amp in self.segments:
            self.segments[amp].set_dark_pix(dark_pix[amp], frac)

    def add_bright_pix(self, bright_pix, nsig=5):
        for amp in self.segments:
            self.segments[amp].add_bright_pix(bright_pix[amp], nsig=nsig)

    def add_traps(self, ndefects, cycles, trap_size):
        traps = OrderedDict()
        for amp in self.segments:
            traps[amp] = self.segments[amp].add_traps(ndefects, cycles,
                                                      trap_size)
        return traps

    def writeto(self, outfile, pars=None, bitpix=-32, obs_time=None):
        ccd_segments = [self.segments[amp] for amp in self.segments]
        output = fitsFile(ccd_segments)
        if pars is not None:
            output[0].header['CCDGAIN'] = pars.system_gain
            output[0].header['BIASLVL'] = pars.bias_level
            output[0].header['CCDNOISE'] = pars.bias_sigma
            output[0].header['RDNOISE'] = pars.read_noise
            output[0].header['DARKCURR'] = pars.dark_current
        for key, value in self.md.items():
            output[0].header[key] = value
        if bitpix > 0:
            my_round = np.round
        else:
            #
            # Delete any BSCALE and BZERO entries, since we are
            # writing image data as floats.
            #
            try:
                del output[amp].header['BSCALE']
                del output[amp].header['BZERO']
            except KeyError:
                pass

            def my_round(x): return x
        for hdu in output[1:-2]:
            hdu.data = np.array(my_round(hdu.data), dtype=self.dtypes[bitpix])
        if obs_time is None:
            # Compute the start of the observation from the current time
            # minus the exposure time.
            obs_time = utcnow(dt=-output[0].header['EXPTIME'])
        output[0].header['DATE-OBS'] = obs_time.isot
        output[0].header['DATE'] = obs_time.isot
        output[0].header.set('MJD-OBS', value=float('%.5f' % obs_time.mjd))
        fitsWriteto(output, outfile, clobber=True, checksum=True)


class SegmentExposure(object):
    def __init__(self, exptime=1, gain=5, ccdtemp=-95, full_well=None,
                 geometry=AmplifierGeometry()):
        self.exptime = exptime
        self.gain = gain
        self.ccdtemp = ccdtemp
        self.full_well = full_well
        self.geometry = geometry
        self.fe55_yield = Fe55Yield(ccdtemp)
        self.image = afwImage.ImageF(geometry.full_segment)
        self.imarr = self.image.Factory(self.image, geometry.imaging).getArray()
        self.ny, self.nx = self.imarr.shape
        self.npix = self.nx*self.ny
        self._sigma = -1

    def add_bias(self, level=1e4, sigma=4):
        """The parameters level and bias are in units of e- and
        converted on output to DN via the system gain."""
        fullarr = self.image.getArray()
        ny, nx = fullarr.shape
        bias_arr = np.array(random.normal(level, sigma, nx*ny),
                            dtype=np.float).reshape(ny, nx)
        fullarr += bias_arr/self.gain

    def add_dark_current(self, level=2e-3):
        """Units of level should be e- per unit time and converted to
        DN on output."""
        dark_arr = self._poisson_imarr(level*self.exptime)/self.gain
        self.imarr += dark_arr

    def expose_flat(self, intensity=0):
        """The parameter intensity is in units of incident photons per
        pixel per unit time."""
        Ne_pred = intensity*self.exptime  # Number of e-/pixel assuming QE=1.
        flat_arr = self._poisson_imarr(Ne_pred)/self.gain
        self.imarr += flat_arr
        if self.full_well is not None:
            indx = np.where(self.imarr > self.full_well/self.gain)
            self.imarr[indx] = self.full_well/self.gain

    def _poisson_imarr(self, Ne):
        return random.poisson(Ne, self.npix).reshape(self.ny, self.nx)

    def sigma(self):
        if self._sigma == -1:
            self._sigma = np.std(self.imarr)
        return self._sigma

    def generate_bright_cols(self, ncols=1):
        bright_cols = np.arange(self.nx)
        random.shuffle(bright_cols)
        return bright_cols[:ncols]

    def set_dark_cols(self, columns, frac):
        for i in columns:
            self.imarr[:, i] *= frac

    def add_bright_cols(self, columns, nsig=5):
        for i in columns:
            self.imarr[:, i] += nsig*self.sigma()

    def generate_bright_pix(self, npix=100):
        bright_pix = np.concatenate((np.ones(npix), np.zeros(self.npix-npix)))
        random.shuffle(bright_pix)
        bright_pix = bright_pix.reshape(self.ny, self.nx)
        return bright_pix

    def set_dark_pix(self, dark_pix, frac):
        self.imarr -= self.imarr*dark_pix*(1. - frac)

    def add_bright_pix(self, bright_pix, nsig=5):
        self.imarr += bright_pix*nsig*self.sigma()

    def add_traps(self, ndefects, cycles, trap_size):
        traps = self.generate_bright_pix(npix=ndefects)
        y, x = np.where(traps != 0)
        xy_pairs = []
        for i, j in zip(x, y):
            # Avoid depleted pixel below first row of imaging area
            if j == 0:
                continue
            self.imarr[j-1][i] -= cycles*(trap_size/self.gain)
            self.imarr[j][i] += cycles*(trap_size/self.gain)
            xy_pairs.append((i, j))
        return xy_pairs

    def add_Fe55_hits(self, nxrays=1000, beta_frac=0.12, sigma=None):
        """
        If sigma is not None, then the e- will be dispersed following a
        2D normal distribution.  sigma is expressed in units of linear
        pixel size.
        """
        Ne_alpha = self.fe55_yield.alpha()[0]
        Ne_beta = self.fe55_yield.beta()[0]
        ny, nx = self.imarr.shape
        if sigma is None:
            # Generated charge per hit is contained within a single pixel.
            for i in range(nxrays):
                x0 = random.randint(nx)
                y0 = random.randint(ny)
                if random.random() < beta_frac:
                    signal = Ne_beta/self.gain
                else:
                    signal = Ne_alpha/self.gain
                self.imarr[y0][x0] += signal
        else:
            # Draw interaction point from full imaging region and e-
            # pixel distribution from 2D Gaussian.
            x0_values = random.uniform(0, nx, nxrays)
            y0_values = random.uniform(0, ny, nxrays)
            for x0, y0 in zip(x0_values, y0_values):
                if random.random() < beta_frac:
                    Ne = Ne_beta
                else:
                    Ne = Ne_alpha
                xvals = random.normal(x0, sigma, size=Ne)
                yvals = random.normal(y0, sigma, size=Ne)
                for xx, yy in zip(xvals, yvals):
                    try:
                        self.imarr[yy][xx] += 1./self.gain
                    except IndexError:
                        pass
        # Round DN/pixel to nearest integer.
        self.imarr = np.round(self.imarr)

    def add_spot_image(self, dn, xref, yref, radius):
        r2 = radius**2
        for x in range(xref-radius, xref+radius, 1):
            for y in range(yref-radius, yref+radius, 1):
                if ((x - xref)**2 + (y - yref)**2 < r2):
                    self.imarr[y][x] += dn

    def add_sys_xtalk_col(self, dn, column):
        self.imarr[:, column] += dn


def fitsFile(ccd_segments):
    headers = fits_headers()
    output = fits.HDUList()
    output.append(fits.PrimaryHDU())
    output[0].header = headers['PRIMARY'].copy()
    output[0].header["EXPTIME"] = ccd_segments[0].exptime
    output[0].header["CCDTEMP"] = ccd_segments[0].ccdtemp
    for amp, segment in zip(imutils.allAmps(), ccd_segments):
        output.append(fits.ImageHDU(data=segment.image.getArray()))
        output[amp].header = headers[headers.keys()[amp]].copy()
        output[amp].header['BZERO'] = 0
        output[amp].name = 'Segment%s' % imutils.channelIds[amp]
        output[amp].header['DETSIZE'] = segment.geometry[amp]['DETSIZE']
        output[amp].header['DATASEC'] = segment.geometry[amp]['DATASEC']
        output[amp].header['DETSEC'] = segment.geometry[amp]['DETSEC']
        output[amp].header['CHANNEL'] = amp
    # Add Test Condition and CCD Operating Condition headers with dummy info.
    output.append(fits.ImageHDU())
    for keyword in headers['TEST_COND']:
        if keyword not in output[-1].header.keys():
            output[-1].header.set(keyword, headers['TEST_COND'][keyword])
    output.append(fits.ImageHDU())
    for keyword in headers['CCD_COND']:
        if keyword not in output[-1].header.keys():
            output[-1].header.set(keyword, headers['CCD_COND'][keyword])
    return output


def writeFits(ccd_segments, outfile, clobber=True):
    output = fitsFile(ccd_segments)
    if clobber:
        try:
            os.remove(outfile)
        except OSError:
            pass
    fitsWriteto(output, outfile, clobber=clobber, checksum=True)
    return outfile


def simulateDark(outfile, dark_curr, exptime=1, hdus=16, verbose=True):
    if verbose:
        print("simulating dark:", outfile)
    segments = []
    for i in range(hdus):
        if verbose:
            print("HDU", i)
        seg = SegmentExposure(exptime=exptime)
        seg.add_bias()
        seg.add_dark_current(dark_curr)
        segments.append(seg)
    writeFits(segments, outfile)


def simulateFlat(outfile, level, gain, dark_curr=1, exptime=1, hdus=16,
                 verbose=True):
    if verbose:
        print("simulating flat:", outfile)
    segments = []
    for i in range(hdus):
        if verbose:
            print("HDU", i)
        seg = SegmentExposure(exptime=exptime, gain=gain)
        seg.add_bias()
        seg.add_dark_current(dark_curr)
        seg.expose_flat(level)
        segments.append(seg)
    writeFits(segments, outfile)


if __name__ == '__main__':
    seg = SegmentExposure()

    seg.add_bias(1e4, 10)
    seg.add_dark_current(1e-2)
    seg.expose_flat(200)
#    cols = seg.add_bright_cols(ncols=1, nsig=5)
#    pix = seg.add_bright_pix(npix=100, nsig=10)
    seg.add_spot_image(2000, 250, 250, 20)

    writeFits((seg,), 'test_image.fits')

    ds9.mtv(seg.image)
