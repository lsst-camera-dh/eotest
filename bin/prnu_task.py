#!/usr/bin/env python
"""
@brief Pixel response non-uniformity.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import lsst.eotest.sensor as sensorTest

parser = sensorTest.TaskParser('Compute pixel response non-uniformity')
parser.add_argument('-f', '--file_pattern', type=str,
                    help='file pattern for PRNU input files')
parser.add_argument('-F', '--prnu_file_list', type=str,
                    help='list of files to be used in PRNU calculations')
parser.add_argument('-c', '--correction_image', type=str,
                    help='correction image for illumination non-uniformity',
                    default=None)
args = parser.parse_args()

task = sensorTest.PrnuTask()
task.config.output_dir = args.output_dir
task.config.verbose = args.verbose

infiles = args.files(args.file_pattern, args.prnu_file_list)
task.run(args.sensor_id, infiles, args.mask_files(), args.system_gains(),
         args.correction_image)
