"""
@brief Task to find traps from pocket-pumped exposure.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import lsst.eotest.image_utils as imutils
from MaskedCCD import MaskedCCD
from EOTestResults import EOTestResults
from Traps import Traps
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
    def run(self, sensor_id, pocket_pumped_file, mask_files, gains,
            cycles=100, threshold=200):
        if self.config.verbose:
            self.log.info("processing %s" % pocket_pumped_file)
        ccd = MaskedCCD(pocket_pumped_file, mask_files=mask_files)
        outfile = os.path.join(self.config.output_dir,
                               '%s_traps.fits' % sensor_id)
        my_traps = Traps(ccd, gains, cycles=cycles)
        my_traps.write(outfile, clobber=True)
        results_file = self.config.eotest_results_file
        if results_file is None:
            results_file = os.path.join(self.config.output_dir,
                                        '%s_eotest_results.fits' % sensor_id)
        results = EOTestResults(results_file)
        if self.config.verbose:
            self.log.info("Amp     Number of traps")
        for amp in imutils.allAmps:
            #
            # Tabulate forward traps (positive polarity) with size >=
            # threshold.
            #
            forward_traps = [item for item in my_traps[amp] 
                             if (item[-1] > 0 and item[-2] >= threshold)]
            num_traps = len(forward_traps)
            results.add_seg_result(amp, 'NUM_TRAPS', num_traps)
            if self.config.verbose:
                self.log.info("%i             %i" % (amp, num_traps))
        results.write(clobber=True)
