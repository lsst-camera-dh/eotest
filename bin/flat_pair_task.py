#!/usr/bin/env python
"""
@brief Compute detector response vs incident flux for pairs of flats from
which linearity and full-well are computed.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import lsst.eotest.sensor as sensorTest

parser = sensorTest.TaskParser('Compute detector response vs incident flux')
parser.add_argument('-f', '--flats', type=str,
                    help='flat pairs file pattern')
parser.add_argument('-F', '--flats_file_list', type=str,
                    help='list of flat pair files')
parser.add_argument('-D', '--detector_response_file', type=str,
                    help='file containing detector response vs incident flux',
                    default=None)
args = parser.parse_args()

task = sensorTest.FlatPairTask()
task.config.output_dir = args.output_dir
task.config.verbose = args.verbose

flat_files = args.files(args.flats, args.flats_file_list)

bias_frame = args.bias_frame('%s_flat_bias_frame.fits' % args.sensor_id)
task.run(args.sensor_id, flat_files, args.mask_files(flat_files[0]),
         args.system_gains(), detrespfile=args.detector_response_file,
         bias_frame=bias_frame)
