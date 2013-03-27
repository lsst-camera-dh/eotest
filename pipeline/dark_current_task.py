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
import image_utils as imUtils
import pipeline.pipeline_utils as pipeUtils

def write_dark_current_maps(images, dark95s, darks, outfile):
    output = pyfits.HDUList()
    output.append(pyfits.PrimaryHDU())
    for i, dark in enumerate(darks):
        output[0].header.update('DARK%02i' % i, os.path.basename(dark))
    for amp in imUtils.allAmps:
        output.append(pyfits.ImageHDU(data=images[amp].getArray()))
        output[0].header.update("DARK95%s" % imUtils.channelIds[amp],
                                dark95s[amp])
        output[amp].name = 'AMP%s' % imUtils.channelIds[amp]
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
        sensor_id = os.environ['SENSOR_ID']
        darks_list = '%s_DARK.txt' % sensor_id
        darks = [x.strip() for x in open(darks_list)]
        outputdir = os.environ['OUTPUTDIR']

    try:
        os.makedirs(outputdir)
    except OSError:
        pass

    gains, sensor = pipeUtils.setup(sys.argv, 5)

    md = afwImage.readMetadata(darks[0], 1)
    exptime = md.get('EXPTIME')

    #
    # Check tempertures
    #
    ccd_temps = [afwImage.readMetadata(x, 1).get('CCDTEMP') for x in darks]
    temp_avg = imUtils.mean(ccd_temps)
    tol = 1.5
    if max(ccd_temps) - temp_avg > tol or temp_avg - min(ccd_temps) > tol:
        raise RuntimeError("Temperature deviations > %s " % tol +
                           "deg C relative to average.")

    def median_image(darks, amp):
        fitsmedian = imUtils.fits_median(darks, imUtils.dm_hdu(amp))
        return imUtils.unbias_and_trim(fitsmedian)
    
    medims = {}
    dark95s = {}
    print "Segment    95 percentile    median"
    for amp in imUtils.allAmps:
        medims[amp] = median_image(darks, amp)
        medims[amp] *= gains[amp]/exptime
        imarr = medims[amp].getArray()
        pixels = imarr.reshape(1, imarr.shape[0]*imarr.shape[1])[0]
        pixels.sort()
        dark95s[amp] = pixels[len(pixels)*0.95]
        sensor.add_seg_result(amp, 'darkCurrent95', dark95s[amp])
        print "%s         %.2e         %.2e" % (imUtils.channelIds[amp],
                                                dark95s[amp],
                                                pixels[len(pixels)/2])

    dark95mean = np.mean(dark95s.values())
    print "CCD: mean 95 percentile value =", dark95mean
    sensor.add_ccd_result('darkCurrent95mean', dark95mean)

    outfile = '%s_dark_current_map.fits' % sensor_id.replace('-', '_')
    outfile = os.path.join(outputdir, outfile)
    write_dark_current_maps(medims, dark95s, darks, outfile)
