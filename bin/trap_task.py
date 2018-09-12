#!/usr/bin/env python
"""
@brief Task to find traps from pocket-pumped exposure.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import glob
import lsst.eotest.sensor as sensorTest

parser = sensorTest.TaskParser('Find Charge Traps')
parser.add_argument('-f', '--pocket_pumped_file', type=str,
                    help='Pocket pumped file')
parser.add_argument('-O', '--output_file', type=str,
                    help='Output text file to contain trap locations and magnitudes')
parser.add_argument('--cycles', type=int, default=100,
                    help='Number of pocket pumping cycles')
parser.add_argument('--threshold', type=int, default=200,
                    help='Trap size threshold (electrons)')
parser.add_argument('--C2_thresh', type=int, default=10,
                    help='C2 correlator detection threshold for trap candidates')
parser.add_argument('--C3_thresh', type=int, default=1,
                    help='C3 detection threshold for trap candidates')
args = parser.parse_args()

task = sensorTest.TrapTask()
task.config.output_dir = args.output_dir
task.config.verbose = args.verbose
task.config.C2_thresh = args.C2_thresh
task.config.C3_thresh = args.C3_thresh

pocket_pumped_file = glob.glob(args.pocket_pumped_file)[0]

task.run(args.sensor_id, pocket_pumped_file,
         args.mask_files(pocket_pumped_file), args.system_gains(),
         cycles=args.cycles, threshold=args.threshold)
