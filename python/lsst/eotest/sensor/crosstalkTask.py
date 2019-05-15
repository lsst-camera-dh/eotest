from __future__ import print_function
from __future__ import absolute_import
import os
import lsst.eotest.image_utils as imutils
from .crosstalk import make_stamp, find_aggressors, crosstalk_model_fit, CrosstalkMatrix
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
                                     int, default=5)
    nsig = pexConfig.Field("Outlier rejection threshold in number of standard deviations of stamp.",
                           float, default=2.0)
    output_dir = pexConfig.Field("Output directory", str, default='.')
    output_file = pexConfig.Field("Output filename", str, default=None)
    verbose = pexConfig.Field("Turn verbosity on", bool, default=True)

class CrosstalkTask(pipeBase.Task):
    """Task to evaluate crosstalk within a single CCD."""
    ConfigClass = CrosstalkConfig
    _DefaultName = "CrosstalkTask"

    @pipeBase.timeMethod
    def run(self, sensor_id, infiles, gains, bias_frame=None, **kwargs):
        #
        # Check kwargs for separate victim sensor information
        #
        infiles2 = None
        sensor_id2 = sensor_id
        gains2 = gains
        if 'infiles2' in kwargs:
            infiles2 = kwargs['infiles2']
            sensor_id2 = kwargs['sensor_id2']
            gains2 = kwargs['gains2']
            bias_frame2 = kwargs.get('bias_frame2', None)
            
        nsig = self.config.nsig
        num_iter = self.config.num_iterations
        threshold = self.config.threshold
        crosstalk_matrix = CrosstalkMatrix(sensor_id, victim_sensor_id=sensor_id2)
        #
        # Perform crosstalk calculation for each aggressor/victim file pair
        #
        for n, infile in enumerate(infiles):
            ccd = MaskedCCD(infile, bias_frame=bias_frame)

            if infiles2 is None:
                victim_ccd = MaskedCCD(infile, bias_frame=bias_frame)
            else:
                victim_ccd = MaskedCCD(infiles2[n], bias_frame=bias_frame2)
            #
            # Aggressor files may have multiple aggressor spots
            # 
            amp_info = find_aggressors(ccd, threshold=threshold)
            for amp, y, x in amp_info:

                row = dict()
                stamp = make_stamp(ccd, amp, y, x)*gains[amp]
                for victim_amp in imutils.allAmps():
                    victim_stamp = make_stamp(victim_ccd, victim_amp, y, x)*gains2[victim_amp]
                    try:
                        row[victim_amp] = crosstalk_model_fit(stamp, victim_stamp,
                                                              num_iter=num_iter, nsig=nsig)
                    except Exception as err:
                        print(err)
                        print("Error during crosstalk calculation.")
                        print("Skipping...")
                crosstalk_matrix.set_row(amp, row)
        #
        # Save crosstalk results to file
        #
        output_dir = self.config.output_dir
        if self.config.output_file is None:
            output_file = os.path.join(output_dir,
                                       '{0}_{1}_crosstalk_results.fits'.format(sensor_id, sensor_id2))
        else:
            output_file = os.path.join(output_dir, self.config.output_file)
        if self.config.verbose:
            self.log.info("Writing crosstalk results file to {0}".format(output_file))
        crosstalk_matrix.write_fits(output_file)
