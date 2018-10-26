#!/usr/bin/env python

"""
@brief Dark current task: compute 95th percentile dark current in
units of e-/sec/pixel.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
from __future__ import print_function
import lsst.eotest.sensor as sensorTest

parser = sensorTest.TaskParser('Compute 95th percentile dark current.')
parser.add_argument('-f', '--dark_files', type=str,
                    help='file pattern for darks')
parser.add_argument('-F', '--dark_file_list', type=str,
                    help='file contain list of dark files')
parser.add_argument('--median_stack', default=None, type=str,
                    help='apply a bias correction')
args = parser.parse_args()

task = sensorTest.DarkCurrentTask()
task.config.output_dir = args.output_dir
task.config.verbose = args.verbose

dark_files = args.files(args.dark_files, args.dark_file_list)
if args.verbose:
    print('processing files:')
    for item in dark_files:
        print('  ', item)

bias_frame = args.bias_frame('%s_dark_bias_frame.fits' % args.sensor_id)
task.run(args.sensor_id, dark_files, args.mask_files(dark_files[0]),
         args.system_gains(), bias_frame=bias_frame)
