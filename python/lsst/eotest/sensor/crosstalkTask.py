"""
@brief Task to produce crosstalk matrix from a set of spot image files.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import lsst.eotest.image_utils as imutils
from lsst.eotest.sensor.crosstalk import make_crosstalk_matrix, CrosstalkMatrix
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

class CrosstalkConfig(pexConfig.Config):
    """Configuration for CrosstalkTask"""
    temp_set_point = pexConfig.Field("Required temperature (C) set point",
                                     float, default=-95.)
    temp_set_point_tol = pexConfig.Field("Required temperature set point tolerance (degrees C)",
                                         float, default=1.)
    output_dir = pexConfig.Field("Output directory", str, default=".")
    verbose = pexConfig.Field("Turn verbosity on", bool, default=True)

class CrosstalkTask(pipeBase.Task):
    """Task to evaluate crosstalk within a single CCD."""
    ConfigClass = CrosstalkConfig
    _DefaultName = "CrosstalkTask"

    @pipeBase.timeMethod
    def run(self, sensor_id, xtalk_files, mask_files, system_xtalk_file=None):
        imutils.check_temperatures(xtalk_files, self.config.temp_set_point_tol,
                                   setpoint=self.config.temp_set_point,
                                   warn_only=True)
        #
        # Test if we have a single (and therefore multi-aggressor) file.
        #
        if len(xtalk_files) == 1:
            xtalk_files = xtalk_files[0]
        xtalk = make_crosstalk_matrix(xtalk_files, mask_files=mask_files)
        if system_xtalk_file is not None:
            system_xtalk_matrix = CrosstalkMatrix(system_xtalk_file)
        xtalk = xtalk - system_xtalk_matrix
        xtalk.write_fits(os.path.join(self.config.output_dir,
                                      '%s_xtalk_matrix.fits' % sensor_id))
