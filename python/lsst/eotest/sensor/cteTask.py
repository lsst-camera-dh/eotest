"""
@brief Compute charge transfer (in)efficiency.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import numpy as np
import astropy.io.fits as fits
from lsst.eotest.fitsTools import fitsWriteto
import lsst.eotest.image_utils as imutils
from AmplifierGeometry import makeAmplifierGeometry
from EOTestResults import EOTestResults
from eperTask import EPERTask

import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

def superflat(files, bias_files=(), outfile='superflat.fits', bitpix=-32,
              bias_subtract=True):
    """
    The superflat is created by bias-offset correcting the input files
    and median-ing them together.
    """
    if bias_files:
        bias_frame = 'mean_bias_frame.fits'
        imutils.fits_mean_file(bias_files, outfile=bias_frame, bitpix=bitpix)

    # Use the first file as a template for the fits output.
    output = fits.open(files[0])
    for amp in imutils.allAmps(files[0]):
        images = afwImage.vectorImageF()
        for infile in files:
            image = afwImage.ImageF(infile, imutils.dm_hdu(amp))
            if bias_subtract:
                if bias_files:
                    bias_image = afwImage.ImageF(bias_frame,
                                                 imutils.dm_hdu(amp))
                else:
                    geom = makeAmplifierGeometry(infile)
                    overscan = geom.serial_overscan
                    bias_image = imutils.bias_image(image,
                                                    overscan=overscan,
                                                    statistic=np.median)
                image -= bias_image
            images.push_back(image)
        median_image = afwMath.statisticsStack(images, afwMath.MEDIAN)
        output[amp].data = median_image.getArray()
        if bitpix is not None:
            imutils.set_bitpix(output[amp], bitpix)
    fitsWriteto(output, outfile, clobber=True)
    return outfile

class CteConfig(pexConfig.Config):
    """Configuration for charge transfer efficiency task"""
    overscans = pexConfig.Field("Number of overscan rows/columns to use",
                                int, default=2)
    output_dir = pexConfig.Field("Output directory", str, default=".")
    eotest_results_file = pexConfig.Field('EO test results filename',
                                          str, default=None)
    verbose = pexConfig.Field("Turn verbosity on", bool, default=True)

class CteTask(pipeBase.Task):
    """Charge transfer efficiency task"""
    ConfigClass = CteConfig
    _DefaultName = "CteTask"

    @pipeBase.timeMethod
    def run(self, sensor_id, superflat_files, bias_files=(), flux_level='high',
            gains=None):
        if flux_level not in ('high', 'low'):
            raise RuntimeError('CteTask: flux_level must be "high" or "low"')
        if self.config.verbose:
            self.log.info("Processing superflat files:")
            for item in superflat_files:
                self.log.info(item)
        #
        # Prepare the co-added superflat file.  Bias subtraction is
        # handled in eperTask.py.
        #
        outfile = '%(sensor_id)s_superflat_%(flux_level)s.fits' % locals()
        superflat_file = superflat(superflat_files, bias_files,
                                   outfile=outfile, bias_subtract=False)
        all_amps = imutils.allAmps(superflat_file)
        #
        # Compute serial CTE.
        #
        s_task = EPERTask()
        s_task.config.direction = 's'
        s_task.config.verbose = self.config.verbose
        s_task.config.cti = True
        scti, bias_ests = s_task.run(superflat_file, all_amps,
                                     self.config.overscans, gains=gains)
        #
        # Compute parallel CTE.
        #
        p_task = EPERTask()
        p_task.config.direction = 'p'
        p_task.config.verbose = self.config.verbose
        p_task.config.cti = True
        pcti, bias_ests = p_task.run(superflat_file, all_amps,
                                     self.config.overscans, gains=gains)
        #
        # Write results to the output file.
        #
        results_file = self.config.eotest_results_file
        if results_file is None:
            results_file = os.path.join(self.config.output_dir,
                                        '%s_eotest_results.fits' % sensor_id)
        results = EOTestResults(results_file, namps=len(all_amps))
        if self.config.verbose:
            self.log.info('CTE %s flux level results' % flux_level)
            self.log.info('amp  parallel_cti  serial_cti')
        for amp in all_amps:
            line = '%i  %s  %s' % (amp, pcti[amp], scti[amp])
            results.add_seg_result(amp, 'CTI_%s_SERIAL' % flux_level.upper(), 
                                   scti[amp].value)
            results.add_seg_result(amp, 'CTI_%s_PARALLEL' % flux_level.upper(),
                                   pcti[amp].value)
            results.add_seg_result(amp,
                                   'CTI_%s_SERIAL_ERROR' % flux_level.upper(), 
                                   scti[amp].error)
            results.add_seg_result(amp, 
                                   'CTI_%s_PARALLEL_ERROR' % flux_level.upper(),
                                   pcti[amp].error)
            if self.config.verbose:
                self.log.info(line)
        results.write(clobber='yes')
