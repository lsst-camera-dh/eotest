"""
@brief Dark current task: compute 95th percentile dark current in
units of e-/sec/pixel.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import sys
import glob
import numpy as np
import pyfits
import lsst.afw.image as afwImage
from image_utils import fits_median, allAmps, dm_hdu, unbias_and_trim, \
     channelIds, mean
from bright_pix import BrightPix
from pipeline.pipeline_utils import setup

def writeFits(images, outfile):
    output = pyfits.HDUList()
    output.append(pyfits.PrimaryHDU())
    for amp in allAmps:
        output.append(pyfits.ImageHDU(data=images[amp].getArray()))
        output[amp].name = 'AMP%s' % channelIds[amp]
    output.writeto(outfile, clobber=True)

def write_darks(darks, outfile):
    output = open(outfile, 'w')
    for dark in darks:
        output.write("%s\n" % dark)
    output.close()

if __name__ == '__main__':
    if len(sys.argv) >= 4:
        darks = glob.glob(sys.argv[1])
        darks.sort()
        darks_list = 'darks.txt'
        write_darks(darks, darks_list)
        sensor_id = sys.argv[2]
        outputdir = sys.argv[3]
    else:
        darks = [x.strip() for x in open(os.environ['DARKS_LIST'])]
        outputdir = os.environ['OUTPUTDIR']
        sensor_id = os.environ['SENSOR_ID']

    try:
        os.makedirs(outputdir)
    except OSError:
        pass

    gains, sensor = setup(sys.argv, 5)

    exptime = afwImage.readMetadata(darks[0], 1).get('EXPTIME')

    #
    # Check tempertures
    #
    ccd_temps = [afwImage.readMetadata(x, 1).get('CCDTEMP') for x in darks]
    temp_avg = mean(ccd_temps)
    tol = 1.5
    if max(ccd_temps) - temp_avg > tol or temp_avg - min(ccd_temps) > tol:
        raise RuntimeError("Temperature deviations > %s " % tol +
                           "deg C relative to average.")

    median_images = {}
    for amp in allAmps:
        median_images[amp] = fits_median(darks, dm_hdu(amp))

    medfile = 'median_dark.fits'
    writeFits(median_images, medfile)

    bright_pix = BrightPix()

    results = bright_pix(medfile, allAmps)

    sensor.add_ccd_result('numBrightPixels', results[0])
    print "Total bright pixels:", results[0]
    print "Segment     # bright pixels"
    for amp, count in zip(allAmps, results[1]):
        sensor.add_seg_result(amp, 'numBrightPixels', count)
        print "%s          %i" % (channelIds[amp], count)
