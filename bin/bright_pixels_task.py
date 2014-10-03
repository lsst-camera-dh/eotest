#!/usr/bin/env python

"""
@brief Bright pixels task: Find pixels and columns in a median image
constructed from an ensemble of darks.  The bright pixel threshold is
specified via the --ethresh option and is in units of -e per pixel per
second.  The threshold for the number of bright pixels that define a
bright column is specified via the --colthresh option.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import lsst.eotest.sensor as sensorTest

parser = sensorTest.TaskParser('Find bright pixels and columns')
parser.add_argument('-f', '--dark_files', type=str,
                    help='file pattern for darks')
parser.add_argument('-F', '--dark_file_list', type=str,
                    help='file containing list of dark files')
parser.add_argument('-e', '--ethresh', default=5, type=int,
                    help='bright pixel threshold in e- per pixel per second')
parser.add_argument('-c', '--colthresh', default=20, type=int,
                    help='bright column threshold in # of bright pixels')
parser.add_argument('-p', '--mask_plane', default='BAD', type=str,
                    help='mask plane to be used for output mask file')
parser.add_argument('-t', '--temp_tol', default=1.5, type=float,
                    help='temperature tolerance for CCDTEMP among dark files')
args = parser.parse_args()

task = sensorTest.BrightPixelsTask()

task.config.ethresh = args.ethresh
task.config.colthresh = args.colthresh
task.config.mask_plane = args.mask_plane
task.config.temp_tol = args.temp_tol
task.config.output_dir = args.output_dir
task.config.eotest_results_file = args.results_file
task.config.verbose = args.verbose

dark_files = args.files(args.dark_files, args.dark_file_list)
if args.verbose:
    print "processing files: "
    for item in dark_files:
        print "  ", item

bias_frame = args.bias_frame('%s_dark_bias_frame.fits' % args.sensor_id)
task.run(args.sensor_id, dark_files, args.mask_files(dark_files[0]), 
         args.system_gains(), bias_frame=bias_frame)
