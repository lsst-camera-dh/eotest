"""
@brief Module to generate sensor-level ratio plots to 
test bias subtraction performance.
"""
from __future__ import print_function
from __future__ import absolute_import

import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import lsst.eotest.image_utils as imutils
from .MaskedCCD import MaskedCCD
from .AmplifierGeometry import makeAmplifierGeometry
import lsst.afw.image as afwImage
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

class SuperbiasConfig(pexConfig.Config):
    """Configuration for SuperbiasTask"""
    output_dir = pexConfig.Field("Output directory", str, default=".")
    eotest_results_file = pexConfig.Field('EO test results filename', 
                                           str, default=None)
    verbose = pexConfig.Field("Turn verbosity on", bool, default=True)

class SuperbiasTask(pipeBase.Task):
    """Task to evaluate bias subtraction via ratio plots."""
    ConfigClass = SuperbiasConfig
    _DefaultName = "SuperbiasTask"

    @pipeBase.timeMethod
    def run(self, sensor_id, sflatL_files, sflatH_files, mask_files, bias_frame):
        
        ## Note: bias_frame is a path to a 16-HDU fits file containing a superbias
        ## for each amplifier

        nl = len(sflatL_files)
        nh = len(sflatH_files)

        ## Make numpy arrays to store the values for the bias and offset-corrected
        ## low/high superflats
        im = MaskedCCD(sflatL_files[0])[1]
        oscan = makeAmplifierGeometry(sflatL_files[0])
        dim = imutils.trim(im, oscan.imaging).getImage().getArray().shape
        lo_num = np.empty((nl, dim[0], dim[1]), dtype=np.float32)
        hi_den = np.empty((nh, dim[0], dim[1]), dtype=np.float32)      
        
        ratio_images = {}
        every_superbias = MaskedCCD(bias_frame, mask_files=mask_files)

        fig, axs = plt.subplots(4,4, figsize=(20,20))
        axs = axs.ravel()

        for amp in imutils.allAmps(sflatL_files[0]):

            superbias = every_superbias[amp]

            ## Calculate the numerator
            for i in range(nl):
                im_lo = MaskedCCD(sflatL_files[i], mask_files=mask_files)[amp]
                im_lo_unbiased = imutils.unbias_and_trim(im=im_lo, overscan=oscan.serial_overscan, 
                                            imaging=oscan.imaging, bias_frame=superbias, 
                                            bias_method='row')
                lo_num[i] = im_lo_unbiased.getImage().getArray()

            ## Calculate the denominator
            for j in range(nh):
                im_hi = MaskedCCD(sflatH_files[j], mask_files=mask_files)[amp]
                im_hi_unbiased = imutils.unbias_and_trim(im=im_hi, overscan=oscan.serial_overscan, 
                                            imaging=oscan.imaging, bias_frame=superbias,
                                            bias_method='row')
                hi_den[j] = im_hi_unbiased.getImage().getArray()

            ## Take the average over all images    
            num = np.mean(lo_num, axis=0) ## Could also try using median stack?
            den = np.mean(hi_den, axis=0)
            ratio = num / den
            ratio_images[amp] = afwImage.ImageF(ratio)
            trimmed_bias_frame = imutils.trim(superbias, oscan.imaging).getImage().getArray()

            ## Limit the range of the histogram to exclude values that fall within the lowest/highest 5% 
            lo_s, hi_s = np.percentile(trimmed_bias_frame, (0.05, 99.95))
            lo_r, hi_r = np.percentile(ratio, (0.05, 99.95))

            ## Plot
            idx = amp - 1
            img = axs[idx].hist2d(trimmed_bias_frame.flatten(), ratio.flatten(), 
                                  bins=([400,200]), range=((lo_s, hi_s), (lo_r, hi_r)),
                                  norm=mpl.colors.LogNorm())
            cbar = fig.colorbar(img[-1], ax=axs[idx])
            axs[idx].set_title('Amp {}'.format(amp))
            if idx%4 == 0:
                axs[idx].set_ylabel('ratio (DN)')

        axs[12].set_xlabel('superbias (DN)')
        axs[13].set_xlabel('superbias (DN)')
        axs[14].set_xlabel('superbias (DN)')
        axs[15].set_xlabel('superbias (DN)')
        plt.suptitle(sensor_id, fontsize=16, y=1.02)
        plt.tight_layout()

        ## Save pdf of all amplifiers
        pdf_file = os.path.join(self.config.output_dir, 
                                '%s_hist2d_plots.pdf' % sensor_id)
        plt.savefig(pdf_file)

        ## File containing the ratio data for each amplifier
        ratio_file = os.path.join(self.config.output_dir, 
                                 '%s_ratio_images.fits' % sensor_id)
        imutils.writeFits(ratio_images, ratio_file, sflatL_files[0], bitpix=-32)

