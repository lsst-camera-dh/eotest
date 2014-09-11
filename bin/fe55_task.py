#!/usr/bin/env python
import os
import lsst.eotest.sensor as sensorTest
import lsst.eotest.image_utils as imutils

parser = sensorTest.TaskParser('PSF and system gain characterization from Fe55 data')
parser.add_argument('-f', '--file_pattern', type=str,
                    help='file pattern for Fe55 input files')
parser.add_argument('-F', '--Fe55_file_list', type=str,
                    help='file name of list of Fe55 files')
parser.add_argument('-c', '--chiprob_min', type=float, default=0.1,
                    help='Mininum chi-square probability for cluster fit')
parser.add_argument('-n', '--nsig', type=float, default=4,
                    help='Footprint threshold in bias noise stdevs')
parser.add_argument('-O', '--outfile', type=str, default=None,
                    help='Output file name (basename only). Computed value if left at default of None: <SENSOR_ID>_psf_results_nsig<nsig>.fits')
parser.add_argument('-x', '--fit_xy', action='store_true', default=False,
                    help='Flag if Gaussian width is to be fit separately in x- and y-directions')
                    
args = parser.parse_args()

task = sensorTest.Fe55Task()

task.config.chiprob_min = args.chiprob_min
task.config.output_dir = args.output_dir
task.config.output_file = args.outfile
task.config.nsig = args.nsig
task.config.verbose = args.verbose
task.config.eotest_results_file = args.results_file
task.config.fit_xy = args.fit_xy

infiles = args.files(args.file_pattern, args.Fe55_file_list)

bias_frame = args.bias_frame('%s_fe55_bias_frame.fits' % args.sensor_id)

task.run(args.sensor_id, infiles, args.mask_files(infiles[0]),
         bias_frame=bias_frame)

os.remove(bias_frame)
