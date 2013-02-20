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
from pipeline.pipeline_utils import setup

def write_dark_current_maps(outfile, images, dark95s, darks):
    output = pyfits.HDUList()
    output.append(pyfits.PrimaryHDU())
    for i, dark in enumerate(darks):
        output[0].header.update('DARK%02i' % i, os.path.basename(dark))
    for amp in allAmps:
        output.append(pyfits.ImageHDU(data=images[amp].getArray()))
        output[0].header.update("DARK95%s" % channelIds[amp], dark95s[amp])
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

    md = afwImage.readMetadata(darks[0], 1)
    exptime = md.get('EXPTIME')

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
    dark95s = {}
    print "Segment    95 percentile dark current"
    for amp in allAmps:
        median_images[amp] = unbias_and_trim(fits_median(darks, dm_hdu(amp)))
        median_images[amp] *= gains[amp]/exptime
        imarr = median_images[amp].getArray()
        pixels = imarr.reshape(1, imarr.shape[0]*imarr.shape[1])[0]
        pixels.sort()
        dark95s[amp] = pixels[len(pixels)*0.95]
        sensor.add_seg_result(amp, 'darkCurrent95', dark95s[amp])
        print "%s         %.3f" % (channelIds[amp], dark95s[amp])

    dark95mean = np.mean(dark95s.values())
    print "CCD: mean 95 percentile value =", dark95mean
    sensor.add_ccd_result('darkCurrent95mean', dark95mean)

    outfile = '%s_dark_current_map.fits' % sensor_id.replace('-', '_')
    outfile = os.path.join(outputdir, outfile)
    write_dark_current_maps(outfile, median_images, dark95s, darks)
