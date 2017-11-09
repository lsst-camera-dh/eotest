"""
@brief PSF characterization from distribution of Gaussian fit
parameters to Fe55 data.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
from MaskedCCD import MaskedCCD
from pipeline.TaskParser import TaskParser
import pylab_plotter as plot
from fe55_psf import PsfGaussFit, afwMath

if __name__ == '__main__':
    parser = TaskParser('PSF characterization')
    parser.add_argument('-f', '--file_pattern', type=str,
                        help='file pattern for Fe55 input files')
    parser.add_argument('-F', '--Fe55_file_list', type=str,
                        help='file name of list of Fe55 files')
    args = parser.parse_args()

    sensor_id = args.sensor_id

    files = args.files(args.file_pattern, args.Fe55_file_list)

    fitter = PsfGaussFit()
    for infile in files:
        print os.path.basename(infile)
        ccd = MaskedCCD(infile, mask_files=args.mask_files())
        for amp in ccd:
            print "   amp", amp
            fitter.process_image(ccd, amp)
    outfile = os.path.join(args.output_dir, '%s_psf_params.fits' % sensor_id)
    fitter.write_results(outfile=outfile)

    sigma, dn, chiprob, amp = fitter.results()

    flags = afwMath.MEDIAN | afwMath.STDEVCLIP

    stats = afwMath.makeStatistics(sigma, flags)
    median = stats.getValue(afwMath.MEDIAN)
    stdev = stats.getValue(afwMath.STDEVCLIP)
    plot.histogram(sigma, xname='Fitted sigma values',
                   xrange=(median-3*stdev, median+3*stdev))
    plot.pylab.savefig(os.path.join(args.output_dir,
                                    '%s_sigma_distribution.png' % sensor_id))
