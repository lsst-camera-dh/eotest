#!/usr/bin/env python
"""
@brief Compute read noise distributions for a sample of images.  Bias
and system readout noise exposures, the latter for determining the
noise contribution from the electronics, must be provided.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import lsst.eotest.sensor as sensorTest

parser = sensorTest.TaskParser('Compute Read Noise')
parser.add_argument('-n', '--noise', type=str, 
                    help='system noise file pattern', default=None)
parser.add_argument('-N', '--noise_file_list', type=str,
                    help='list of system noise files', default=None)
parser.add_argument('-x', '--dx', default=100, type=int,
                    help='subregion size in pixels along x-direction')
parser.add_argument('-y', '--dy', default=100, type=int,
                    help='subregion size in pixels along y-direction')
parser.add_argument('-S', '--nsamp', default=1000, type=int,
                    help='number of subregions to sample')
parser.add_argument('--use_overscan', default=False, action='store_true',
                    help='Use serial overscan region for noise estimates.')
args = parser.parse_args()

task = sensorTest.ReadNoiseTask()
task.config.dx = args.dx
task.config.dy = args.dy
task.config.nsamp = args.nsamp
task.config.output_dir = args.output_dir
task.config.eotest_results_file = args.results_file

bias_files = args.files(args.bias_frame_pattern, args.bias_frame_list)
if args.noise is None and args.noise_file_list is None:
    system_noise_files = None
else:
    system_noise_files = args.files(args.noise, args.noise_file_list)

task.run(args.sensor_id, bias_files, args.system_gains(),
         system_noise_files=system_noise_files,
         mask_files=args.mask_files(bias_files[0]),
         use_overscan=args.use_overscan)
