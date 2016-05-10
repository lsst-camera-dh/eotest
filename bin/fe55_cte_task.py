#!/usr/bin/env python

"""
Characterize CTE using Fe55 cluster asymmetries
"""
import lsst.eotest.sensor as sensorTest

parser = sensorTest.TaskParser('Characterize CTE using Fe55 cluster asymmetries')
parser.add_argument('-f', '--fe55_file_pattern', type=str,
                    help='Fe55 file pattern')
parser.add_argument('-F', '--fe55_file_list', type=str,
                    help='list of Fe55 files')
parser.add_argument('--readout_direction', type=str, default='serial',
                    help='readout direction to use: serial or parallel')

args = parser.parse_args()

task = sensorTest.Fe55CteTask()
task.config.output_dir = args.output_dir
task.config.direction = args.readout_direction
task.config.verbose = args.verbose

fe55_files = args.files(args.fe55_file_pattern, args.fe55_file_list)
task.run(args.sensor_id, fe55_files, args.mask_files(fe55_files[0]))
