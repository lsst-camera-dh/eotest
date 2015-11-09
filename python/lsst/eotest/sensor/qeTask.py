"""
@brief Compute QE curves from the wavelength scan dataset.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import QE
import lsst.eotest.image_utils as imutils
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

class QeConfig(pexConfig.Config):
    """Configuration for QE measurement task"""
    temp_set_point = pexConfig.Field("Required temperature (C) set point",
                                     float, default=-95.)
    temp_set_point_tol = pexConfig.Field("Required temperature set point tolerance (degrees C)",
                                         float, default=1.)
    output_dir = pexConfig.Field("Output directory", str, default=".")
    verbose = pexConfig.Field("Turn verbosity on", bool, default=True)

class QeTask(pipeBase.Task):
    """Task to compute QE curves from wavelength scan dataset"""
    ConfigClass = QeConfig
    _DefaultName = "QeTask"

    @pipeBase.timeMethod
    def run(self, sensor_id, qe_files, pd_ratio_file, mask_files, gains,
            medians_file=None, vendor_data=False, correction_image=None):
        imutils.check_temperatures(qe_files, self.config.temp_set_point_tol,
                                   setpoint=self.config.temp_set_point,
                                   warn_only=True)
        qe_data = QE.QE_Data(verbose=self.config.verbose, logger=self.log)
        
        if medians_file is None:
            medians_file = os.path.join(self.config.output_dir,
                                        '%s_QE_medians.txt' % sensor_id)
            qe_data.calculate_medians(qe_files, medians_file,
                                      mask_files=mask_files, clobber=True,
                                      correction_image=correction_image)
            
        qe_data.read_medians(medians_file)
        
        if vendor_data:
            qe_data.incidentPower_e2v()
        else:
            qe_data.incidentPower(pd_ratio_file)

        qe_data.calculate_QE(gains)

        fits_outfile = os.path.join(self.config.output_dir,
                                    '%s_QE.fits' % sensor_id)
        qe_data.write_fits_tables(fits_outfile)
