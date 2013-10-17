#!/usr/bin/env python
"""
@brief Task to find traps from pocket-pumped exposure.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import lsst.eotest.sensor as sensorTest

parser = sensorTest.TaskParser('Find Charge Traps')
parser.add_argument('-f', '--pocket_pumped_file', type=str,
                    help='Pocket pumped file')
parser.add_argument('-O', '--output_file', type=str,
                    help='Output text file to contain trap locations and magnitudes')
args = parser.parse_args()

task = sensorTest.TrapTask()
task.config.output_dir = args.output_dir
task.config.verbose = args.verbose

task.run(args.sensor_id, args.pocket_pumped_file, args.mask_files(),
         args.system_gains())
