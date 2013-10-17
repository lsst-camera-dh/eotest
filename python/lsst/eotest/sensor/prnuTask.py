"""
@brief Pixel response non-uniformity.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
from MaskedCCD import Metadata
from prnu import prnu
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

class PrnuConfig(pexConfig.Config):
    """Configuration for pixel response non-uniformity task"""
    output_dir = pexConfig.Field("Output directory", str, default=".")
    verbose = pexConfig.Field("Turn verbosity on", bool, default=True)

class PrnuTask(pipeBase.Task):
    """Task for computing pixel response non-uniformity"""
    ConfigClass = PrnuConfig
    _DefaultName = "PrnuTask"

    @pipeBase.timeMethod
    def run(self, sensor_id, prnu_files, mask_files, gains, correction_image):
        outfile = os.path.join(self.config.output_dir,
                               '%s_prnu_values.txt' % sensor_id)
        output = open(outfile, 'w')
        line = "wavelength (nm)   pixel_stdev   pixel_median"
        output.write(line + '\n')
        if self.config.verbose:
            self.log.info(line)
        for infile in prnu_files:
            md = Metadata(infile, 1)
            wl = md.get('MONOWL')
            if int(wl) in (350, 450, 500, 620, 750, 870, 1000):
                pix_stdev, pix_median = prnu(infile, mask_files, gains,
                                             correction_image=correction_image)
                line = "%6.1f  %12.4e  %12.4e" % (wl, pix_stdev, pix_median)
                output.write(line + '\n')
                if self.config.verbose:
                    self.log.info(line)
        output.close()
