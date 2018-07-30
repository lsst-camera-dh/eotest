import os
import lsst.eotest.image_utils as imutils
from lsst.eotest.sensor.crosstalk_new import make_crosstalk_matrix, CrosstalkMatrix
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

def bias_subtracted_image(image, bias_image, overscan, fit_order=1,
                          statistic=np.median):
    # Make deep copies of image and bias image so that we can modify them.
    im_out = image.Factory(image)
    bias_sub = bias_image.Factory(bias_image)
    # Subtract overscans.
    bias_sub -= imutils.bias_image(bias_image, overscan, fit_order=fit_order,
                                   statistic=statistic)
    im_out -= imutils.bias_image(image, overscan, fit_order=fit_order,
                                 statistic=statistic)
    # Subtract remaining strucutured bias.
    im_out -= bias_sub
    return im_out

def stackedxtalk(xtalk_files, bias_frame=None, outfile='stackedxtalk.fits', 
                 bias_subtract=True):
    
    overscan = makeAmplifierGeometry(files[0]).serial_overscan
    output = fits.open(files[0])
    for amp in imutils.allAmps(files[0]):
        images = afwImage.vectorImageF()
        for infile in files:
            image = afwImage.ImageF(infile, imutils.dm_hdu(amp))
            if bias_subtract:
                if bias_frame:
                    bias_image = afwImage.ImageF(bias_frame, 
                                                 imutils.dm_hdu(amp))
                    image = bias_subtracted_image(image, bias_image, overscan)
                else:
                    image -= imutils.bias_image(image, overscan, 
                                                statistics=np.median)

            images.push_back(image)
        median_image = afwMath.statisticsStack(images, afwMath.MEDIAN)
        output[amp].data = median_image.getArray()
        if bitpix is not None:
            imutils.set_bitpix(output[amp], bitpix)
    fitsWriteto(output, outfile, clobber=True)
    return outfile

class CrosstalkConfig(pexConfig.Config):
    
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
    def run(self, sensor_id, xtalk_files, mask_files, bias_frame=None,
            gains=None, system_xtalk_file=None):
        imutils.check_temperatures(xtalk_file, self.config.temp_set_point_tol,
                                   setpoint=self.config.temp_set_point,
                                   warn_only=True)

        positions = set([x.split('_')[-3] for x in xtalk_files])
        stackedxtalk_files = []
        for pos in positions:
            position_files = [x for x in xtalk_files if pos in x]
            outfile = '%(sensor_id)s_stackedxtalk_%(pos)s.fits' % locals()
            stackedxtalk_file = stackedxtalk(position_files, bias_frame=bias_frame,
                                             outfile=outfile)
            stackedxtalk_files.append(stackedxtalk_file)

        xtalk = make_crosstalk_matrix(stackedxtalk_files, mask_files=mask_files)
        if system_xtalk_file is not None:
            system_xtalk_matrix = CrosstalkMatrix(systtem_xtalk_file)
        xtalk = xtalk - system_xtalk_matrix
        xtalk.write_fits(os.path.join(self.config.output_dir, 
                                      '%s_xtalk_matrix.fits' % sensor_id)
                                         
