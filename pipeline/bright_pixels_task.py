"""
@brief Bright pixels task: Find pixels and columns in a median image
constructed from an ensemble of darks.  The brightness threshold is
specified in nsig of the noise above the mean.

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
import simulation.sim_tools as sim_tools

from bright_pix import BrightPix

def writeFits(images, outfile):
    output = pyfits.HDUList()
    output.append(pyfits.PrimaryHDU())
    for amp in imUtils.allAmps:
        output.append(pyfits.ImageHDU(data=images[amp].getArray()))
        output[amp].name = 'AMP%s' % imUtils.channelIds[amp]
        output[amp].header.update('DETSIZE', imUtils.detsize)
        output[amp].header.update('DETSEC', imUtils.detsec(amp))
    output.writeto(outfile, clobber=True)

def write_darks_list(darks, outfile):
    output = open(outfile, 'w')
    for dark in darks:
        output.write("%s\n" % dark)
    output.close()

if __name__ == '__main__':
    if len(sys.argv) >= 4:
        darks = glob.glob(sys.argv[1])
        darks.sort()
        darks_list = 'darks.txt'
        write_darks_list(darks, darks_list)
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

    exptime = afwImage.readMetadata(darks[0], 1).get('EXPTIME')

    #
    # Check tempertures
    #
    ccd_temps = [afwImage.readMetadata(x, 1).get('CCDTEMP') for x in darks]
    temp_avg = imUtils.mean(ccd_temps)
    tol = 1.5
    if max(ccd_temps) - temp_avg > tol or temp_avg - min(ccd_temps) > tol:
        raise RuntimeError("Temperature deviations > %s " % tol +
                           "deg C relative to average.")

    median_images = {}
    for amp in imUtils.allAmps:
        median_images[amp] = imUtils.fits_median(darks, imUtils.dm_hdu(amp))

    medfile = os.path.join(outputdir, '%s_median_dark_bp.fits' % sensor_id)
    writeFits(median_images, medfile)

    try:
        nsig = int(os.environ['NSIG'])
    except KeyError:
        nsig = 5

    bright_pix = BrightPix(nsig=nsig)

    results = bright_pix(medfile, imUtils.allAmps)

    sensor.add_ccd_result('numBrightPixels', results[0])
    print "Total bright pixels:", results[0]
    print "Segment     # bright pixels"
    for amp, count in zip(imUtils.allAmps, results[1]):
        sensor.add_seg_result(amp, 'numBrightPixels', count)
        print "%s          %i" % (imUtils.channelIds[amp], count)

    #
    # Create images of bright pixels and columns.
    #
    segments = []
    for amp in imUtils.allAmps:
        seg = sim_tools.SegmentExposure()
        for ix, iy in results[2][amp-1]:
            seg.imarr[iy, ix] = 1
        segments.append(seg)
        for ix in results[3][amp-1]:
            seg.imarr[:, ix] = 1
    outfile = os.path.join(outputdir, '%s_bright_pixel_map.fits' % sensor_id)
    sim_tools.writeFits(segments, outfile)
