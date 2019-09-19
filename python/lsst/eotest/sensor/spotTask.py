"""
@brief Source identification and characterization of
spot images.
"""
from __future__ import print_function
from __future__ import absolute_import
import os
import numpy as np
from astropy.io import fits

import lsst.eotest.image_utils as imutils
import lsst.afw.image as afwImage
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.pipe.tasks
from lsst.pipe.tasks.characterizeImage import CharacterizeImageTask, CharacterizeImageConfig
import lsst.meas.extensions.shapeHSM

from .MaskedCCD import MaskedCCD
from .AmplifierGeometry import parse_geom_kwd
    
def make_ccd_mosaic(infile, bias_frame=None, dark_frame=None, gains=None):
    """Create a full CCD image mosaic.
    
    Combine the 16 amplifier arrays into a single array, performing
    bias/offset subtraction and gain correction.
    
    Args:
        infile: Image FITs filename.
        bias_frame: Bias frame FITs filename.
        dark_frame: Dark frame FITs filename.
        gains: Dictionary mapping amplifier number to gain.
        
    Returns:
        afwImage object containing mosaiced image.
    """
    
    ccd = MaskedCCD(infile, bias_frame=bias_frame)

    ## Additional images for calibration
    if dark_frame is not None:
        dark = MaskedCCD(dark_frame, bias_frame=bias_frame)

    ## Get amp geometry information
    foo = fits.open(infile)
    datasec = parse_geom_kwd(foo[1].header['DATASEC'])
    nx_segments = 8
    ny_segments = 2
    nx = nx_segments*(datasec['xmax'] - datasec['xmin'] + 1)
    ny = ny_segments*(datasec['ymax'] - datasec['ymin'] + 1)

    mosaic = np.zeros((ny, nx), dtype=np.float32)

    for ypos in range(ny_segments):
        for xpos in range(nx_segments):
            amp = ypos*nx_segments + xpos + 1

            detsec = parse_geom_kwd(foo[amp].header['DETSEC'])
            xmin = nx - max(detsec['xmin'], detsec['xmax'])
            xmax = nx - min(detsec['xmin'], detsec['xmax']) + 1
            ymin = ny - max(detsec['ymin'], detsec['ymax'])
            ymax = ny - min(detsec['ymin'], detsec['ymax']) + 1

            ## Dark correction
            if dark_frame is not None:
                exptime = ccd.md.get('EXPTIME')
                dark_exptime = dark.md.get('EXPTIME')

                imarr = ccd.unbiased_and_trimmed_image(amp).getImage().getArray()
                dark_imarr = dark.unbiased_and_trimmed_image(amp).getImage().getArray()
                subarr = imarr - dark_imarr*exptime/dark_exptime
            else:
                subarr = ccd.unbiased_and_trimmed_image(amp).getImage().getArray()
                    
            ## Flip array orientation (if applicable)
            if detsec['xmax'] > detsec['xmin']: # flip in x-direction
                subarr = subarr[:, ::-1]
            if detsec['ymax'] > detsec['ymin']: # flip in y-direction
                subarr = subarr[::-1, :]

            ## Gain correction
            if gains is not None:
                subarr *= gains[amp]

            ## Assign to final array
            mosaic[ymin:ymax, xmin:xmax] = subarr

    image = afwImage.ImageF(mosaic)
    return image

class SpotConfig(pexConfig.Config):
    """Configuration for Spot analysis task"""
    minpixels = pexConfig.Field("Minimum number of pixels above detection threshold", 
                                int, default=10)
    nsig = pexConfig.Field("Source footprint threshold in number of standard deviations.", 
                           float, default=10)
    bgbinsize = pexConfig.Field("Bin size for background estimation.", int, default=10)
    output_dir = pexConfig.Field("Output directory", str, default='.')
    output_file = pexConfig.Field("Output filename", str, default=None)
    verbose = pexConfig.Field("Turn verbosity on", bool, default=True)

class SpotTask(pipeBase.Task):
    """Task to estimate spot moments from spot projector data."""

    ConfigClass = SpotConfig
    _DefaultName = "SpotTask"

    @pipeBase.timeMethod
    def run(self, sensor_id, infile, gains, bias_frame=None, flat_frame=None, dark_frame=None):
        #
        # Process a CCD image mosaic
        #
        if self.config.verbose:
            self.log.info("processing {0}".format(infile))
        image = make_ccd_mosaic(infile, bias_frame=bias_frame, 
                                dark_frame=dark_frame, gains=gains)
        if flat_frame is not None:
            flat = make_ccd_mosaic(flat_frame, bias_frame=bias_frame,
                                   dark_frame=dark_frame, gains=gains)
            flatarr = flat.getArray()
            imarr = image.getArray()
            image = afwImage.ImageF(imarr*np.median(flatarr)/flatarr)
        exposure = afwImage.ExposureF(image.getBBox())
        exposure.setImage(image)
        #
        # Set up characterize task configuration
        #
        nsig = self.config.nsig
        bgbinsize = self.config.bgbinsize
        minpixels = self.config.minpixels
        charConfig = CharacterizeImageConfig()
        charConfig.doMeasurePsf = False
        charConfig.doApCorr = False
        charConfig.repair.doCosmicRay = False
        charConfig.detection.minPixels = minpixels
        charConfig.detection.background.binSize = bgbinsize
        charConfig.detection.thresholdType = "stdev"
        charConfig.detection.thresholdValue = nsig
        hsm_plugins = set(["ext_shapeHSM_HsmShapeBj",
                           "ext_shapeHSM_HsmShapeLinear",
                           "ext_shapeHSM_HsmShapeKsb",
                           "ext_shapeHSM_HsmShapeRegauss",
                           "ext_shapeHSM_HsmSourceMoments",
                           "ext_shapeHSM_HsmPsfMoments"])  
        charConfig.measurement.plugins.names |= hsm_plugins
        charTask = CharacterizeImageTask(config=charConfig)
        result = charTask.run(exposure)

        src = result.sourceCat
        if self.config.verbose:
            self.log.info("Detected {0} objects".format(len(src)))
        #
        # Save catalog results to file
        #
        output_dir = self.config.output_dir
        if self.config.output_file is None:
            output_file = os.path.join(output_dir,
                                       '{0}_source_catalog.cat'.format(sensor_id))
        else:
            output_file = os.path.join(output_dir, self.config.output_file)
        if self.config.verbose:
            self.log.info("Writing spot results file to {0}".format(output_file))
        src.writeFits(output_file)
