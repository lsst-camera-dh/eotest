"""
@brief Persistence task.  Compute the deferred charge as a function of time
in darks taken after a saturated flat has been taken.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import numpy as np
import astropy.io.fits as fits
from lsst.eotest.fitsTools import fitsTableFactory, fitsWriteto
import lsst.eotest.image_utils as imutils
from MaskedCCD import MaskedCCD
from EOTestResults import EOTestResults
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
import lsst.pex.exceptions.wrappers as pexExceptWrap
import lsst.pipe.base as pipeBase
import datetime
import astropy.time

def file_timestamp(infile):
    timestamp = infile[-len('YYYYMMDDHHMMSS.fits'):-len('.fits')]
    year, month, day, hours, minutes, seconds = \
        [int(timestamp[:4])] + [int(timestamp[4+2*i:4+2*(i+1)])
                                for i in range(5)]
    t_readout = astropy.time.Time(datetime.datetime(year, month, day, hours, 
                                                    minutes, seconds),
                                  format='datetime', scale='utc')
    return t_readout

def readout_time(infile):
    """
    Return readout time as an astropy.time.Time object.
    It is computed from the PHDU keyword values as MJD-OBS + EXPTIME
    """
    md = afwImage.readMetadata(infile, 1)
    mjd_obs = astropy.time.Time(md.get('MJD-OBS'), format='mjd',
                                scale='utc')
    t_readout = astropy.time.Time(mjd_obs.datetime + 
                                  datetime.timedelta(seconds=md.get('EXPTIME')),
                                  scale='utc')
    return t_readout

class PersistenceConfig(pexConfig.Config):
    """Configuration for PersistenceTask"""
    region_size = pexConfig.Field("Linear size of (square) region-of-interest in pixels",
                                  int, default=200)
    region_x_offset = pexConfig.Field("Pixel offset in x-direction from center of imaging region of a segment", int, default=0)
    region_y_offset = pexConfig.Field("Pixel offset in y-direction from center of imaging region of a segment", int, default=0)
    temp_set_point = pexConfig.Field("Required temperature (C) set point",
                                     float, default=-95.)
    temp_set_point_tol = pexConfig.Field("Required temperature set point tolerance (degrees C)",
                                         float, default=1.)
    output_dir = pexConfig.Field("Output directory", str, default=".")
    eotest_results_file = pexConfig.Field('EO test results filename',
                                          str, default=None)
    verbose = pexConfig.Field("Turn verbosity on", bool, default=True)

class PersistenceTask(pipeBase.Task):
    """Task to evaluate dark current quantiles."""
    ConfigClass = PersistenceConfig
    _DefaultName = "PersistenceTask"

    @pipeBase.timeMethod
    def run(self, sensor_id, pre_flat_darks, flat, post_flat_darks, 
            mask_files, gains):
        darks = list(pre_flat_darks) + list(post_flat_darks)
        imutils.check_temperatures(darks, self.config.temp_set_point_tol,
                                   setpoint=self.config.temp_set_point,
                                   warn_only=True)
        # Check that pre-flat dark frames all have the same exposure time
        md = imutils.Metadata(pre_flat_darks[0], 1)
        exptime = md.get('EXPTIME')
        for item in pre_flat_darks[1:]:
            md = imutils.Metadata(item, 1)
            if exptime != md.get('EXPTIME'):
                raise RuntimeError("Exposure times of pre-flat darks differ.")

        # Make a median image of the preflat darks
        median_images = {}
        for amp in imutils.allAmps(darks[0]):
            median_images[amp] = imutils.fits_median(pre_flat_darks,
                                                     imutils.dm_hdu(amp))
        medfile = os.path.join(self.config.output_dir,
                               '%s_persistence_dark_median.fits' % sensor_id)
        imutils.writeFits(median_images, medfile, darks[0])
        ccd = MaskedCCD(medfile, mask_files=mask_files)

        # Define the sub-region for assessing the deferred charge.
        # This is the same bounding box for all segments, so use amp=1.
        image = ccd.unbiased_and_trimmed_image(1)
        xllc = ((image.getWidth() - self.config.region_size)/2. 
                - self.config.region_x_offset)
        yllc = ((image.getHeight() - self.config.region_size)/2. 
                - self.config.region_y_offset)
        imaging_reg = afwGeom.Box2I(afwGeom.Point2I(int(xllc), int(yllc)),
                                    afwGeom.Extent2I(self.config.region_size,
                                                     self.config.region_size))
        overscan = ccd.amp_geom.serial_overscan
        # Compute reference dark current for each segment.
        dc_ref = {}
        for amp in ccd:
            mi = imutils.unbias_and_trim(ccd[amp], overscan, imaging_reg)
            dc_ref[amp] = afwMath.makeStatistics(mi, afwMath.MEDIAN,
                                                 ccd.stat_ctrl).getValue()
            dc_ref[amp] *= gains[amp]/exptime
        
        # Extract reference time for computing the time dependence
        # of the deferred charge as the observation time + exposure time
        # from the saturated flat.
        tref = readout_time(flat)

        # Loop over post-flat darks, compute median e-/pixel in
        # subregion, subtract dc_ref*exptime, persist, and report the
        # deferred charge vs time (using MJD-OBS + EXPTIME) for each amp.
        deferred_charges = []
        times = []
        for dark in post_flat_darks:
            ccd = MaskedCCD(dark, mask_files=mask_files)
            dt = readout_time(dark) - tref
            times.append(dt.sec)
            exptime = ccd.md.get('EXPTIME')
            charge = {}
            for amp in ccd:
                mi = imutils.unbias_and_trim(ccd[amp], overscan, imaging_reg)
                estimators = afwMath.MEDIAN | afwMath.STDEV
                stats = afwMath.makeStatistics(mi, estimators, ccd.stat_ctrl)
                value = (stats.getValue(afwMath.MEDIAN)*gains[amp] 
                         - dc_ref[amp]*exptime)
                stdev = (stats.getValue(afwMath.STDEV)*gains[amp] 
                         - dc_ref[amp]*exptime)
                charge[amp] = (value, stdev)
            deferred_charges.append(charge)

        if self.config.verbose:
            for amp in ccd:
                self.log.info("amp: %i" % amp)
                for i, time in enumerate(times):
                    self.log.info("%.1f  %e  %e" % 
                                  (time, deferred_charges[i][amp][0], 
                                   deferred_charges[i][amp][1]))

        outfile = os.path.join(self.config.output_dir,
                               '%s_persistence.fits' % sensor_id)
        self.write(times, deferred_charges, outfile, clobber=True)

    @pipeBase.timeMethod
    def write(self, times, deferred_charges, outfile, clobber=True):
        colnames = ['TIME']
        columns = [times]
        units = ['s']
        for amp in deferred_charges[0]:
            colnames.extend(['MEDIAN%02i' % amp, 'STDEV%02i' % amp])
            columns.extend([np.array([deferred_charges[i][amp][0]
                                      for i in range(len(times))]),
                            np.array([deferred_charges[i][amp][1]
                                      for i in range(len(times))])])
            units.extend(['e-/pixel', 'rms e-/pixel'])
        formats = ['E']*len(columns)

        HDUList = fits.HDUList()
        HDUList.append(fits.PrimaryHDU())
        HDUList.append(fitsTableFactory([fits.Column(name=colname,
                                                     format=format,
                                                     unit=unit,
                                                     array=column)
                                           for colname, format, unit, column in
                                           zip(colnames, formats, units, columns)]))
        HDUList[-1].name = 'IMAGE_PERSISTENCE_CURVES'
        fitsWriteto(HDUList, outfile, clobber=clobber)
