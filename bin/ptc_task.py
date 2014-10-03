#!/usr/bin/env python
"""
@brief Compute photon transfer curve for flat pair dataset.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import lsst.eotest.sensor as sensorTest

parser = sensorTest.TaskParser('Compute photon transfer curve')
parser.add_argument('-f', '--flats', type=str,
                    help='flat pairs file pattern')
parser.add_argument('-F', '--flats_file_list', type=str,
                    help='list of flat pair files')
parser.add_argument('-R', '--rebinning_factor', type=int,
                    help='rebinning factor', default=1)
args = parser.parse_args()

task = sensorTest.PtcTask()
task.config.output_dir = args.output_dir
task.config.verbose = args.verbose

flat_files = args.files(args.flats, args.flats_file_list)

bias_frame = args.bias_frame('%s_flat_bias_frame.fits' % args.sensor_id)
task.run(args.sensor_id, flat_files, args.mask_files(flat_files[0]),
         args.system_gains(), binsize=args.rebinning_factor,
         bias_frame=bias_frame)
         
