#!/usr/bin/env python

"""
@brief Compute charge transfer (in)efficiency.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import glob
import pyfits
import lsst.test_scripts.image_utils as imutils
import lsst.test_scripts.sensor as sensorTest
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath

def superflat(files, outfile='superflat.fits'):
    """
    The superflat is created by bias-offset correcting the input files
    and median-ing them together.
    """
    # Use the first file as a template for the pyfits output.
    output = pyfits.open(files[0])
    for amp in imutils.allAmps:
        images = afwImage.vectorImageF()
        for infile in files:
            decorated_image = afwImage.DecoratedImageF(infile,
                                                       imutils.dm_hdu(amp))
            sr = sensorTest.SegmentRegions(decorated_image)
            image = decorated_image.getImage()
            image -= imutils.bias_image(image, overscan=sr.serial_overscan)
            images.push_back(image)
        median_image = afwMath.statisticsStack(images, afwMath.MEDIAN)
        output[amp].data = median_image.getArray()
    output.writeto(outfile, clobber=True)
    return outfile

if __name__ == '__main__':
    parser = sensorTest.TaskParser('Compute charge transfer efficiency in parallel and serial directions using extended pixel edge response technique')
    parser.add_argument('-f', '--superflat_pattern', type=str,
                        help='superflat dataset file pattern')
    parser.add_argument('-F', '--superflat_file_list', type=str,
                        help='list of superflat files')
    args = parser.parse_args()
    sensor = args.sensor()
    sensor_id = args.sensor_id
    outfile = os.path.join(args.output_dir, '%s_cti_values.txt' % sensor_id)

    files = args.files(args.superflat_pattern, args.superflat_file_list)
    superflat_file = superflat(files)

    overscans = 3
    #
    # Compute serial CTE
    #
    s_task = sensorTest.EPERTask()
    s_task.config.direction = 's'
    s_task.config.verbose = False
    s_task.config.cti = True
    scti = s_task.run(superflat_file, imutils.allAmps, overscans)
    #
    # Compute parallel CTE
    #
    p_task = sensorTest.EPERTask()
    p_task.config.direction = 'p'
    p_task.config.verbose = False
    p_task.config.cti = True
    pcti = p_task.run(superflat_file, imutils.allAmps, overscans)
    #
    # Write results to the db and text file
    #
    sensor.add_ccd_result('ctiParallelMean', imutils.mean(pcti.values()))
    sensor.add_ccd_result('ctiSerialMean', imutils.mean(scti.values()))
    output = open(outfile, 'w')
    output.write('amp  parallel_cti  serial_cti\n')
    for amp in imutils.allAmps:
        sensor.add_seg_result(amp, 'ctiParallel', pcti[amp])
        sensor.add_seg_result(amp, 'ctiSerial', scti[amp])
        output.write('%s  %12.4e  %12.4e\n' % (imutils.channelIds[amp],
                                               pcti[amp], scti[amp]))
    output.close()
    #
    # Clean up
    #
    os.remove(superflat_file)
