#!/usr/bin/env python

"""
@brief Compute charge transfer (in)efficiency.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import lsst.eotest.sensor as sensorTest

parser = sensorTest.TaskParser('Compute charge transfer efficiency in parallel and serial directions using extended pixel edge response technique')
parser.add_argument('-f', '--superflat_pattern', type=str,
                    help='superflat dataset file pattern')
parser.add_argument('-F', '--superflat_file_list', type=str,
                    help='list of superflat files')
parser.add_argument('-n', '--overscans', type=int, default=3,
                    help='number of overscan rows/columns to use')

args = parser.parse_args()

task = sensorTest.CteTask()
task.config.overscans = args.overscans
task.config.output_dir = args.output_dir
task.config.verbose = args.verbose

superflat_files = args.files(args.superflat_pattern, args.superflat_file_list)

task.run(args.sensor_id, superflat_files, args.mask_files())
