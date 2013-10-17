"""
@brief Compute charge transfer (in)efficiency.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import pyfits

import lsst.eotest.image_utils as imutils
from MaskedCCD import SegmentRegions
from eperTask import EPERTask

import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

def superflat(files, outfile='superflat.fits'):
    """
    The superflat is created by bias-offset correcting the input files
    and median-ing them together.
    """
    # Use the first file as a template for the pyfits output.
    output = pyfits.open(files[0])
    for amp in imutils.allAmps:
        images = afwImage.vectorImageF()
        for infile in files:
            decorated_image = afwImage.DecoratedImageF(infile,
                                                       imutils.dm_hdu(amp))
            sr = SegmentRegions(decorated_image)
            image = decorated_image.getImage()
            image -= imutils.bias_image(image, overscan=sr.serial_overscan)
            images.push_back(image)
        median_image = afwMath.statisticsStack(images, afwMath.MEDIAN)
        output[amp].data = median_image.getArray()
    output.writeto(outfile, clobber=True)
    return outfile

class CteConfig(pexConfig.Config):
    """Configuration for charge transfer efficiency task"""
    overscans = pexConfig.Field("Number of overscan rows/columns to use",
                                int, default=3)
    output_dir = pexConfig.Field("Output directory", str, default=".")
    verbose = pexConfig.Field("Turn verbosity on", bool, default=True)

class CteTask(pipeBase.Task):
    """Charge transfer efficiency task"""
    ConfigClass = CteConfig
    _DefaultName = "CteTask"

    @pipeBase.timeMethod
    def run(self, sensor_id, superflat_files, mask_files):
        if self.config.verbose:
            self.log.info("Processing superflat files: %s" % superflat_files)
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
                          self.config.overscans)
        #
        # Compute parallel CTE.
        #
        p_task = EPERTask()
        p_task.config.direction = 'p'
        p_task.config.verbose = False
        p_task.config.cti = True
        pcti = p_task.run(superflat_file, imutils.allAmps,
                          self.config.overscans)
        #
        # Write results to the output file.
        #
        outfile = os.path.join(self.config.output_dir,
                               '%s_cti_values.txt' % sensor_id)
        output = open(outfile, 'w')
        output.write('amp  parallel_cti  serial_cti\n')
        if self.config.verbose:
            self.log.info('amp  parallel_cti  serial_cti')
        for amp in imutils.allAmps:
            line = '%s  %12.4e  %12.4e' % (imutils.channelIds[amp],
                                           pcti[amp], scti[amp])
            output.write(line + '\n')
            if self.config.verbose:
                self.log.info(line)
        output.close()
        #
        # Clean up
        #
        os.remove(superflat_file)
