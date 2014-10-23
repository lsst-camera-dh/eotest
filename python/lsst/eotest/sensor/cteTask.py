"""
@brief Compute charge transfer (in)efficiency.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import pyfits
from lsst.eotest.pyfitsTools import pyfitsWriteto
import lsst.eotest.image_utils as imutils
from AmplifierGeometry import makeAmplifierGeometry
from EOTestResults import EOTestResults
from eperTask import EPERTask

import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

def superflat(files, outfile='superflat.fits', bitpix=None):
    """
    The superflat is created by bias-offset correcting the input files
    and median-ing them together.
    """
    # Use the first file as a template for the pyfits output.
    output = pyfits.open(files[0])
    for amp in imutils.allAmps:
        images = afwImage.vectorImageF()
        for infile in files:
            image = afwImage.ImageF(infile, imutils.dm_hdu(amp))
            geom = makeAmplifierGeometry(infile)
            image -= imutils.bias_image(image, overscan=geom.serial_overscan)
            images.push_back(image)
        median_image = afwMath.statisticsStack(images, afwMath.MEDIAN)
        output[amp].data = median_image.getArray()
        if bitpix is not None:
            imutils.set_bitpix(output[amp], bitpix)
    pyfitsWriteto(output, outfile, clobber=True)
    return outfile

class CteConfig(pexConfig.Config):
    """Configuration for charge transfer efficiency task"""
    overscans = pexConfig.Field("Number of overscan rows/columns to use",
                                int, default=3)
    output_dir = pexConfig.Field("Output directory", str, default=".")
    eotest_results_file = pexConfig.Field('EO test results filename',
                                          str, default=None)
    verbose = pexConfig.Field("Turn verbosity on", bool, default=True)

class CteTask(pipeBase.Task):
    """Charge transfer efficiency task"""
    ConfigClass = CteConfig
    _DefaultName = "CteTask"

    @pipeBase.timeMethod
    def run(self, sensor_id, superflat_files, mask_files):
        if self.config.verbose:
            self.log.info("Processing superflat files:")
            for item in superflat_files:
                self.log.info(item)
        #
        # Prepare the co-added superflat file.
        #
        superflat_file = superflat(superflat_files)
        #
        # Compute serial CTE.
        #
        s_task = EPERTask()
        s_task.config.direction = 's'
        s_task.config.verbose = False
        s_task.config.cti = True
        scti = s_task.run(superflat_file, imutils.allAmps,
                          self.config.overscans, mask_files)
        #
        # Compute parallel CTE.
        #
        p_task = EPERTask()
        p_task.config.direction = 'p'
        p_task.config.verbose = False
        p_task.config.cti = True
        pcti = p_task.run(superflat_file, imutils.allAmps,
                          self.config.overscans, mask_files)
        #
        # Write results to the output file.
        #
        results_file = self.config.eotest_results_file
        if results_file is None:
            results_file = os.path.join(self.config.output_dir,
                                        '%s_eotest_results.fits' % sensor_id)
        results = EOTestResults(results_file)
        if self.config.verbose:
            self.log.info('amp  parallel_cti  serial_cti')
        for amp in imutils.allAmps:
            line = '%s  %12.4e  %12.4e' % (imutils.channelIds[amp],
                                           pcti[amp], scti[amp])
            results.add_seg_result(amp, 'CTI_SERIAL', scti[amp])
            results.add_seg_result(amp, 'CTI_PARALLEL', pcti[amp])
            if self.config.verbose:
                self.log.info(line)
        results.write(clobber='yes')
        #
        # Clean up
        #
        os.remove(superflat_file)
