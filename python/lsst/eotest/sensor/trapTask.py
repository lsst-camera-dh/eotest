"""
@brief Task to find traps from pocket-pumped exposure.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
from __future__ import absolute_import, print_function
import os
import lsst.eotest.image_utils as imutils
from .MaskedCCD import MaskedCCD
from .EOTestResults import EOTestResults
from .Traps import Traps
from .generate_mask import generate_mask
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase


class TrapConfig(pexConfig.Config):
    """Configuration for TrapTask"""
    output_dir = pexConfig.Field("Output directory", str, default=".")
    outfile = pexConfig.Field("Output file (base)name", str, default=None)
    C2_thresh = pexConfig.Field("C2 threshold", float, default=10.)
    C3_thresh = pexConfig.Field("C3 threshold", float, default=1.)
    nx = pexConfig.Field("Local background width (pixels)", int, default=10)
    ny = pexConfig.Field("Local background height (pixels)", int, default=10)
    edge_rolloff = pexConfig.Field("Edge rolloff width (pixels)", int,
                                   default=10)
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
        if self.config.outfile is None:
            outfile = os.path.join(self.config.output_dir,
                                   '%s_traps.fits' % sensor_id)
        else:
            outfile = os.path.join(self.config.output_dir,
                                   self.config.outfile)
        my_traps = Traps(ccd, gains, cycles=cycles,
                         C2_thresh=self.config.C2_thresh,
                         C3_thresh=self.config.C3_thresh,
                         nx=self.config.nx, ny=self.config.ny,
                         edge_rolloff=self.config.edge_rolloff)
        my_traps.write(outfile, overwrite=True)
        results_file = self.config.eotest_results_file
        if results_file is None:
            results_file = os.path.join(self.config.output_dir,
                                        '%s_eotest_results.fits' % sensor_id)
        results = EOTestResults(results_file, namps=len(ccd))
        if self.config.verbose:
            self.log.info("Amp     Number of traps")
        columns = {}
        for amp in ccd:
            #
            # Tabulate forward traps (A0 < 0) with size >= threshold.
            #
            forward_traps = [item for item in my_traps[amp]
                             if (item[-2] < 0 and item[-3] >= threshold)]
            columns[amp] = [item[0] for item in forward_traps]
            num_traps = len(forward_traps)
            results.add_seg_result(amp, 'NUM_TRAPS', num_traps)
            if self.config.verbose:
                self.log.info("%i             %i" % (amp, num_traps))
        mask_file = '%s_traps_mask.fits' % sensor_id
        generate_mask(pocket_pumped_file, mask_file, 'TRAPS', columns=columns)
        results.write(clobber=True)
