"""
@brief Task to find traps from pocket-pumped exposure.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import lsst.eotest.image_utils as imutils
from MaskedCCD import MaskedCCD
from traps import traps
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

class TrapConfig(pexConfig.Config):
    """Configuration for TrapTask"""
    output_dir = pexConfig.Field("Output directory", str, default=".")
    eotest_results_file = pexConfig.Field("EO test results filename",
                                          str, default=None)
    verbose = pexConfig.Field("Turn verbosity on", bool, default=True)

class TrapTask(pipeBase.Task):
    """Configuration for task to find traps from pocket-pumped exposure"""
    ConfigClass = TrapConfig
    _DefaultName = "TrapTask"
    @pipeBase.timeMethod
    def run(self, sensor_id, pocket_pumped_file, mask_files, gains):
        if self.config.verbose:
            self.log.info("processing on %s" % pocket_pumped_file)
        ccd = MaskedCCD(pocket_pumped_file, mask_files=mask_files)
        outfile = os.path.join(self.config.output_dir,
                               '%s_traps.txt' % sensor_id)
        my_traps = traps(ccd, gains, outfile=outfile)
        results_file = self.config.eotest_results_file
        if results_file is None:
            resuts_file = '%s_eotest_results.fits' % sensor_id
        results = EOTestResults(results_file)
        for amp in imutils.allAmps:
            results.add_seg_result(amp, 'NUM_TRAPS', len(my_traps[amp]))
        results.write(clobber=True)
