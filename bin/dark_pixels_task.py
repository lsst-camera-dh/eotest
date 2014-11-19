#!/usr/bin/env python

"""
@brief Dark pixels task: Find pixels and columns in a median image
constructed from an ensemble of superflats.  The dark pixel threshold
is specified via the --thresh option and is in units of fraction of
the segment mean.  The threshold for the number of dark pixels that
define a dark column is specified via the --colthresh option.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import lsst.eotest.sensor as sensorTest

parser = sensorTest.TaskParser('Find dark pixels and columns')
parser.add_argument('-f', '--sflat_files', type=str,
                    help='file pattern for superflats')
parser.add_argument('-F', '--sflat_file_list', type=str,
                    help='file containing list of superflat files')
parser.add_argument('-t', '--thresh', default=0.8, type=float,
                    help='fractional dark pixel threshold in units of median pixel value')
parser.add_argument('-c', '--colthresh', default=20, type=int,
                    help='dark column threshold in # of dark pixels')
parser.add_argument('-p', '--mask_plane', default='BAD', type=str,
                    help='mask plane to be used for output mask file')
args = parser.parse_args()

task = sensorTest.DarkPixelsTask()

task.config.thresh = args.thresh
task.config.colthresh = args.colthresh
task.config.mask_plane = args.mask_plane
task.config.output_dir = args.output_dir
task.config.eotest_results_file = args.results_file
task.config.verbose = args.verbose

sflat_files = args.files(args.sflat_files, args.sflat_file_list)
if args.verbose:
    print "processing files: "
    for item in sflat_files:
        print "  ", item

#bias_frame = args.bias_frame('%s_dark_bias_frame.fits' % args.sensor_id)
bias_frame = None
task.run(args.sensor_id, sflat_files, args.mask_files(sflat_files[0]),
         bias_frame=bias_frame)
