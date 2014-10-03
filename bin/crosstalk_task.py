#!/usr/bin/env python
"""
@brief Task to produce crosstalk matrix from a set of spot image files.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import pylab
import lsst.eotest.sensor as sensorTest

parser = sensorTest.TaskParser('Compute crosstalk from a set of spot images')
parser.add_argument('-f', '--xtalk_files', type=str,
                    help='file pattern for crosstalk files')
parser.add_argument('-F', '--xtalk_file_list', type=str,
                    help='list of crosstalk files')

args = parser.parse_args()

task = sensorTest.CrosstalkTask()
task.config.output_dir = args.output_dir
task.config.verbose = args.verbose

xtalk_files = args.files(args.xtalk_files, args.xtalk_file_list)

task.run(args.sensor_id, xtalk_files, args.mask_files(xtalk_files[0]))
