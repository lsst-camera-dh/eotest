#!/usr/bin/env python

"""
@brief PSF characterization and system gain from distribution of
Gaussian fit parameters to Fe55 data.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import numpy as np
import pyfits
import lsst.afw.math as afwMath
import lsst.test_scripts.image_utils as imutils
import lsst.test_scripts.sensor as sensorTest
import lsst.test_scripts.sensor.pylab_plotter as plot

if __name__ == '__main__':
    parser = sensorTest.TaskParser('PSF and system gain characterization from Fe55 data')
    parser.add_argument('-f', '--file_pattern', type=str,
                        help='file pattern for Fe55 input files')
    parser.add_argument('-F', '--Fe55_file_list', type=str,
                        help='file name of list of Fe55 files')
    parser.add_argument('-O', '--output_gain_file', type=str,
                        help='Output FITS file to contain gain values')
    parser.add_argument('-c', '--chiprob_min', type=float, default=0.1,
                        help='Mininum chi-square probability for cluster fit')
    args = parser.parse_args()

    sensor = args.sensor()
    sensor_id = args.sensor_id
    
    files = args.files(args.file_pattern, args.Fe55_file_list)
    
    fitter = sensorTest.PsfGaussFit()
    for infile in files:
        print os.path.basename(infile)
        ccd = sensorTest.MaskedCCD(infile, mask_files=args.mask_files())
        for amp in ccd:
            print "   amp", amp
            fitter.process_image(ccd, amp)
    outfile = os.path.join(args.output_dir, '%s_psf_params.fits' % sensor_id)
    fitter.write_results(outfile=outfile)
    #
    # Plot sigma values for the entire CCD.
    #
    sigma, dn, dn_fp, chiprob, amp = fitter.results()
    flags = afwMath.MEDIAN | afwMath.STDEVCLIP

    stats = afwMath.makeStatistics(sigma, flags)
    median = stats.getValue(afwMath.MEDIAN)
    stdev = stats.getValue(afwMath.STDEVCLIP)
    plot.histogram(sigma, xname='Fitted sigma values',
                   xrange=(median-3*stdev, median+3*stdev))
    plot.pylab.savefig(os.path.join(args.output_dir,
                                    '%s_sigma_distribution.png' % sensor_id))

    #
    # Fit the DN distributions to obtain the system gain per amp
    #
    gains = {}
    results = pyfits.open(outfile)
    for amp in imutils.allAmps:
        dn = np.array(results[amp].data.field('DN'), dtype=np.float)
        chiprob = results[amp].data.field('CHIPROB')
        indx = np.where(chiprob > args.chiprob_min)
        gains[amp] = sensorTest.fe55_gain_fitter(dn[indx], make_plot=False)

    #
    # Write gain results to db table and output file.
    #
    if args.output_gain_file is not None:
        outfile = os.path.join(args.output_dir, args.output_gain_file)
    else:
        outfile = os.path.join(args.output_dir, "%s_gains.fits" % sensor_id)
                           
    output = pyfits.HDUList()
    output.append(pyfits.PrimaryHDU())
    
    gain_median = imutils.median(gains.values())
    sensor.add_ccd_result('gainMedian', gain_median)
    print "Median gain among segments:", gain_median
    print "Segment    gain"
    for amp in imutils.allAmps:
        sensor.add_seg_result(amp, 'gain', gains[amp])
        print "%s         %.4f" % (imutils.channelIds[amp], gains[amp])
        output[0].header.update("GAIN%s" % imutils.channelIds[amp],
                                gains[amp])
    output.writeto(outfile, clobber=True, checksum=True)
