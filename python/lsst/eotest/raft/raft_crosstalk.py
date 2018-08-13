import argparse
import os
import glob
import numpy as np
import itertools

#from lsst.eotest.sensor.crosstalkTask import stackedxtalk, bias_subtracted_image
#from lsst.eotest.sensor.crosstalk import make_crosstalk_matrix

## Needed just for testing
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.eotest.image_utils as imutils
from lsst.eotest.sensor.AmplifierGeometry import makeAmplifierGeometry
import astropy.io.fits as fits
from lsst.eotest.fitsTools import fitsWriteto

## Needed for testing
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

## Needed for testing
def stackedxtalk(files, bias_frames=None, outfile='stackedxtalk.fits', 
                 bitpix=-32, bias_subtract=True, gains=None):
    """Make a median stacked xtalk image."""

    if gains is None:
        gains = dict([(amp, 1) for amp in imutils.allAmps(files[0])])
    
    overscan = makeAmplifierGeometry(files[0]).serial_overscan
    output = fits.open(files[0])
    for amp in imutils.allAmps(files[0]):
        images = afwImage.vectorImageF()
        for infile in files:
            image = afwImage.ImageF(infile, imutils.dm_hdu(amp))
            ## Perform bias/overscan subtraction
            if bias_subtract:
                ## Calculate median bias frame
                if bias_frames:
                    bias_images = afwImage.vectorImageF()
                    for bias_frame in bias_frames:
                        bias_image = afwImage.ImageF(bias_frame,
                                                     imutils.dm_hdu(amp))
                        bias_images.push_back(image)
                    median_bias = afwMath.statisticsStack(bias_images, afwMath.MEDIAN)
                    image = bias_subtracted_image(image, bias_image, overscan)
                ## Overscan correct
                else:
                    image -= imutils.bias_image(image, overscan, 
                                                statistics=np.median)

            images.push_back(image)
        median_image = afwMath.statisticsStack(images, afwMath.MEDIAN)
        output[amp].data = median_image.getArray()*gains[amp]
        if bitpix is not None:
            imutils.set_bitpix(output[amp], bitpix)
    fitsWriteto(output, outfile, clobber=True)
    return outfile

def aggregate_xtalk_results():
    raise NotImplementedError

def which_ccd(x, y):
    """Determine focal CCD from projector x/y information."""

    if x > -21.0 and x < 21.0:
        i = 1
    elif x > 21.0:
        i = 2
    elif x < -21.0:
        i = 0

    if y > -21.0 and y < 21.0:
        j = 1
    elif y > 21.0:
        j = 0
    elif y < -21.0:
        j = 2

    ccd = 'S{0}{1}'.format(i, j)
    return ccd

class CrosstalkHandler():

    def __init__(self, output_dir, seqno_list=None):

        self.xtalk_dict = dict()
        self.output_dir = output_dir
        self.seqno_list = seqno_list
            
    def ingest_ccd_files(self, sensor_id, image_files, bias_files, 
                         mask_files, dark_files=None):
        """Calibrate crosstalk images for each CCD."""

        aggressor_seqno = []
        calibrated_file_dict = dict()

        if self.seqno_list is None:
            self.seqno_list = list(set([f.split('_')[-3] for f in image_files]))
        for seqno in self.seqno_list:
            imfiles = [f for f in image_files if f.split('_')[-3] == seqno]
            bfiles = [f for f in bias_files if f.split('_')[-3] == seqno]
            if dark_files is not None:
                dfiles = [f for f in dark_files if f.split('_')[-3] == seqno]
            else:
                dfiles = None
            outfile = os.path.join(self.output_dir, 
                                   '{0}_stacked_xtalk_{1}.fits'.format(sensor_id, seqno))
            calibrated_file = stackedxtalk(imfiles, bfiles, outfile)
            projector_x = float(imfiles[0].split('_')[-5])
            projector_y = float(imfiles[0].split('_')[-4])
            if which_ccd(projector_x, projector_y) == sensor_id:
                aggressor_seqno.append(seqno)
            calibrated_file_dict[seqno] = calibrated_file

        self.xtalk_dict[sensor_id] = (aggressor_seqno, calibrated_file_dict,
                                      mask_files)

    def measure_raft_xtalk(self):
        """Measure crosstalk for all sensor combinations."""

        ccd_list = sorted(self.xtalk_dict.keys())
        namps = len(ccd_list)*16
        raft_xtalk_result = CrosstalkMatrix(namps=namps)
        
        for i, (aggressor_sensor_id, victim_sensor_id) in enumerate(itertools.product(ccd_list, ccd_list)):

            print "{0} x {1}".format(aggressor_sensor_id, victim_sensor_id)

            aggressor_seqno = self.xtalk_dict[aggressor_sensor_id][0]
            aggressor_files = [self.xtalk_dict[aggressor_sensor_id][1][seqno]
                               for seqno in aggressor_seqno]
            aggressor_masks = self.xtalk_dict[aggressor_sensor_id][2]
            victim_files = [self.xtalk_dict[victim_sensor_id][1][seqno]
                            for seqno in aggressor_seqno]
            victim_masks = self.xtalk_dict[victim_sensor_id][2]

            result = make_crosstalk_matrix(aggressor_files, aggressor_masks, 
                                           victim_file_list=victim_files,
                                           victim_mask_files=victim_masks)
            
            raft_xtalk_result.matrix[(i/9)*16:(i/9)*16+16,
                                     (i%9)*16:(i%9)*16+16] = result.matrix

        raft_xtalk_result.write_fits(os.path.join(self.output_dir,
                                      'raft_xtalk_matrix.fits'))
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('directory')
    args = parser.parse_args()

    main_dir = args.directory

    ccd_names = ['S00', 'S01', 'S02', 'S10', 'S11', 'S12', 'S20', 'S21', 'S22']

    mask_files = ()
    xtalk_handler = CrosstalkHandler('/nfs/slac/g/ki/ki19/lsst/snyder18/LSST',
                                     ['{0:03d}'.format(i) for i in range(36)])
    for ccd in ccd_names:
        
        all_files = glob.glob(os.path.join(main_dir, ccd, '*.fits'))
        image_files = [f for f in all_files if '_xy_stage_' not in f]
        bias_files = [f for f in all_files if 'bias' in f]
        xtalk_handler.ingest_ccd_files(ccd, image_files, bias_files, mask_files)

    xtalk_handler.measure_raft_xtalk()
