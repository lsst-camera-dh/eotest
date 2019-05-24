"""
@brief Compute charge transfer (in)efficiency.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
from __future__ import absolute_import
import os
import astropy.io.fits as fits
from lsst.eotest.fitsTools import fitsWriteto
import lsst.eotest.image_utils as imutils
from .AmplifierGeometry import makeAmplifierGeometry
from .EOTestResults import EOTestResults
from .eperTask import EPERTask
import lsst.afw
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

def bias_subtracted_image(image, bias_image, overscan, bias_method='row'):
    # Make deep copies of image and bias image so that we can modify them.
    im_out = image.Factory(image)
    bias_sub = bias_image.Factory(bias_image)
    # Subtract overscans.
    bias_sub -= imutils.bias_image(im=bias_image, overscan=overscan, bias_method=bias_method)
    im_out -= imutils.bias_image(im=image, overscan=overscan, bias_method=bias_method)
    # Subtract remaining strucutured bias.
    im_out -= bias_sub
    return im_out


def superflat(files, bias_frame=None, outfile='superflat.fits', bitpix=None,
              bias_subtract=True, bias_method='row'):
    """
    The superflat is created by bias-offset correcting the input files
    and median-ing them together.
    """
    # Get overscan region.
    overscan = makeAmplifierGeometry(files[0]).serial_overscan
    output_images = dict()
    for amp in imutils.allAmps(files[0]):
        images = []
        for infile in files:
            image = afwImage.ImageF(infile, imutils.dm_hdu(amp))
            if bias_subtract:
                if bias_frame:
                    bias_image = afwImage.ImageF(bias_frame,
                                                 imutils.dm_hdu(amp))
                    image = bias_subtracted_image(image, bias_image, overscan,
                                                  bias_method)
                else:
                    image -= imutils.bias_image(im=image, overscan=overscan,
                                                bias_method=bias_method)
            images.append(image)
        if lsst.afw.__version__.startswith('12.0'):
            images = afwImage.vectorImageF(images)
        output_images[amp] = afwMath.statisticsStack(images, afwMath.MEDIAN)
    imutils.writeFits(output_images, outfile, files[0])
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
    def run(self, sensor_id, superflat_files, bias_frame=None,
            flux_level='high', gains=None, mask_files=()):
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
        superflat_file = superflat(superflat_files, bias_frame=bias_frame,
                                   outfile=outfile, bias_subtract=False)
        nframes = len(superflat_files)
        all_amps = imutils.allAmps(superflat_file)
        #
        # Compute serial CTE.
        #
        s_task = EPERTask()
        s_task.config.direction = 's'
        s_task.config.verbose = self.config.verbose
        s_task.config.cti = True
        scti, bias_ests = s_task.run(superflat_file, nframes, all_amps,
                                     self.config.overscans, gains=gains,
                                     mask_files=mask_files)
        #
        # Compute parallel CTE.
        #
        p_task = EPERTask()
        p_task.config.direction = 'p'
        p_task.config.verbose = self.config.verbose
        p_task.config.cti = True
        pcti, bias_ests = p_task.run(superflat_file, nframes, all_amps,
                                     self.config.overscans, gains=gains,
                                     mask_files=mask_files)
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
        results.write(clobber=True)
