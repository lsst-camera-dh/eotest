from __future__ import print_function
from __future__ import absolute_import
import os
import lsst.eotest.image_utils as imutils
from .crosstalk import stamp, find_aggressors, crosstalk_model_fit, CrosstalkMatrix
from .MaskedCCD import MaskedCCD

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

class CrosstalkConfig(pexConfig.Config):
    
    temp_set_point = pexConfig.Field("Required temperature (C) set point",
                                     float, default=-95.)
    temp_set_point_tol = pexConfig.Field("Required temperature set point tolerance (degrees C)",
                                         float, default=1.)
    threshold = pexConfig.Field("Flux threshold for aggressor spots.",
                                float, default=100000.0)
    num_iterations = pexConfig.Field("Number of least-square iterations to perform.",
                                     int, default=10)
    nsig = pexConfig.Field("Outlier rejection threshold in number of standard deviations of stamp.",
                           float, default=3.0)
    output_dir = pexConfig.Field("Output directory", str, default='.')
    output_file = pexConfig.Field("Output filename", str, default=None)
    verbose = pexConfig.Field("Turn verbosity on", bool, default=True)

class CrosstalkTask(pipeBase.Task):
    """Task to evaluate crosstalk within a single CCD."""
    ConfigClass = CrosstalkConfig
    _DefaultName = "CrosstalkTask"

    @pipeBase.timeMethod
    def run(self, sensor_id, aggressor_files, bias_frame=None,
            gains=None, victim_files=None, victim_bias_frame=None):

        # Maybe check later?
        if victim_files is not None and len(victim_files) != len(aggressor_files):
            raise ValueError("Number of aggressor and victim files must match!")

        crosstalk_matrix = CrosstalkMatrix()

        nsig = self.config.nsig
        num_iter = self.config.num_iterations
        threshold = self.config.threshold
        #
        # Perform crosstalk calculation for each aggressor/victim file pair
        #
        for n, aggressor_file in enumerate(aggressor_files):
            aggressor_ccd = MaskedCCD(aggressor_file, bias_frame=bias_frame)

            if victim_files is None:
                victim_ccd = MaskedCCD(aggressor_file, bias_frame=bias_frame)
            else:
                victim_ccd = MaskedCCD(victim_files[n], bias_frame=victim_bias_frame)
            #
            # Aggressor files may have multiple aggressor spots
            # 
            aggressor_amp_info = find_aggressors(aggressor_ccd, threshold=threshold)
            for aggressor_amp, y, x in aggressor_amp_info:

                row = dict()
                aggressor_stamp = stamp(aggressor_ccd, aggressor_amp, y, x)
                for amp in imutils.allAmps():
                    victim_stamp = stamp(victim_ccd, amp, y, x)
                    try:
                        row[amp] = crosstalk_model_fit(aggressor_stamp, victim_stamp,
                                                       num_iter=num_iter, nsig=nsig)
                    except Exception as err:
                        print(err)
                        print("Error during crosstalk calculation.")
                        print("Skipping...")
                crosstalk_matrix.set_row(aggressor_amp, row)
        #
        # Save crosstalk results to file
        #
        output_dir = self.config.output_dir
        if self.config.output_file is None:
            output_file = os.path.join(output_dir,
                                       '{0}_crosstalk_results.fits'.format(sensor_id))
        else:
            output_file = os.path.join(output_dir, self.config.output_file)
        if self.config.verbose:
            self.log.info("Writing crosstalk results file to {0}".format(output_file))
        crosstalk_matrix.write_fits(output_file)
