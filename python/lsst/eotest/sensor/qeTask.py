"""
@brief Compute QE curves from the wavelength scan dataset.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import QE
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

class QeConfig(pexConfig.Config):
    """Configuration for QE measurement task"""
    output_dir = pexConfig.Field("Output directory", str, default=".")
    verbose = pexConfig.Field("Turn verbosity on", bool, default=True)

class QeTask(pipeBase.Task):
    """Task to compute QE curves from wavelength scan dataset"""
    ConfigClass = QeConfig
    _DefaultName = "QeTask"

    @pipeBase.timeMethod
    def run(self, sensor_id, qe_files, ccd_cal_file, sph_cal_file,
            wlscan_file, mask_files, gains):
        qe_data = QE.QE_Data(verbose=self.config.verbose, logger=self.log)
        
        medians_file = os.path.join(self.config.output_dir,
                                    '%s_QE.txt' % sensor_id)
        qe_data.calculate_medians(qe_files, medians_file,
                                  mask_files=mask_files, clobber=True)
        qe_data.read_medians(medians_file)
        
        qe_data.calculate_QE(ccd_cal_file, sph_cal_file, wlscan_file, gains)

        fits_outfile = os.path.join(self.config.output_dir,
                                    '%s_QE.fits' % sensor_id)
        qe_data.write_fits_tables(fits_outfile)
