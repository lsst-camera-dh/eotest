#!/usr/bin/env python
"""
@brief Compute QE curves from the wavelength scan dataset.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import lsst.eotest.sensor as sensorTest

parser = sensorTest.TaskParser('Compute QE curves')
parser.add_argument('-f', '--qe_files', type=str,
                    help='wavelength scan file pattern')
parser.add_argument('-F', '--qe_file_list', type=str,
                    help='list of wavelength scan files')
parser.add_argument('-c', '--ccd_cal_file', type=str,
                    help='calibration file for photodiode at CCD location')
parser.add_argument('-i', '--int_sph_cal_file', type=str,
                    help='calibration file for photodiode at integrating sphere')
parser.add_argument('-w', '--wavelength_scan_file', type=str,
                    help='uncorrected wavelength scan file')
args = parser.parse_args()

task = sensorTest.QeTask()
task.config.output_dir = args.output_dir
task.config.verbose = args.verbose

qe_files = args.files(args.qe_files, args.qe_file_list)

task.run(args.sensor_id, qe_files, args.ccd_cal_file, args.int_sph_cal_file,
         args.wavelength_scan_file, args.mask_files(), args.system_gains())
