#!/usr/bin/env python
"""
@brief Compute read noise distributions for a sample of images.  Bias
and system readout noise exposures, the latter for determining the
noise contribution from the electronics, must be provided.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import lsst.eotest.sensor as sensorTest

parser = sensorTest.TaskParser('Compute Read Noise')
parser.add_argument('-b', '--bias', type=str,
                    help='bias file pattern')
parser.add_argument('-B', '--bias_file_list', type=str,
                    help='list of bias files')
parser.add_argument('-n', '--noise', type=str, 
                    help='system noise file pattern')
parser.add_argument('-N', '--noise_file_list', type=str,
                    help='list of system noise files')
parser.add_argument('-x', '--dx', default=100, type=int,
                    help='subregion size in pixels along x-direction')
parser.add_argument('-y', '--dy', default=100, type=int,
                    help='subregion size in pixels along y-direction')
parser.add_argument('-S', '--nsamp', default=1000, type=int,
                    help='number of subregions to sample')
args = parser.parse_args()

task = sensorTest.ReadNoiseTask()
task.config.dx = args.dx
task.config.dy = args.dy
task.config.nsamp = args.nsamp

bias_files = args.files(args.bias, args.bias_file_list)
system_noise_files = args.files(args.noise, args.noise_file_list)

task.run(args.sensor_id, bias_files, system_noise_files, args.mask_files(),
         args.system_gains())
