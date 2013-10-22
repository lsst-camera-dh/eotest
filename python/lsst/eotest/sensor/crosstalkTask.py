"""
@brief Task to produce crosstalk matrix from a set of spot image files.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
from lsst.eotest.sensor.crosstalk import make_crosstalk_matrix
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

class CrosstalkConfig(pexConfig.Config):
    """Configuration for CrosstalkTask"""
    output_dir = pexConfig.Field("Output directory", str, default=".")
    verbose = pexConfig.Field("Turn verbosity on", bool, default=True)

class CrosstalkTask(pipeBase.Task):
    """Task to evaluate crosstalk within a single CCD."""
    ConfigClass = CrosstalkConfig
    _DefaultName = "CrosstalkTask"

    @pipeBase.timeMethod
    def run(self, sensor_id, xtalk_files, mask_files):
        #
        # Test if we have a single (and therefore multi-aggressor) file.
        #
        if len(xtalk_files) == 1:
            xtalk_files = xtalk_files[0]
        xtalk = make_crosstalk_matrix(xtalk_files, mask_files=mask_files)
        xtalk.write(os.path.join(self.config.output_dir,
                                 '%s_xtalk_matrix.txt' % sensor_id))
