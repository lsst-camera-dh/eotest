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
parser.add_argument('-p', '--pd_ratio_file', type=str,
                    help='photodiode ratio file')
parser.add_argument('-M', '--medians_file', type=str,
                    help='file of median pixel values from wavelength scan dataset', default=None)
parser.add_argument('-c', '--correction_image', type=str,
                    help='correction image for illumination non-uniformity',
                    default=None)
parser.add_argument('--vendor_data', default=False, action='store_true',
                    help='Have data from vendor')
args = parser.parse_args()

task = sensorTest.QeTask()
task.config.output_dir = args.output_dir
task.config.verbose = args.verbose

qe_files = args.files(args.qe_files, args.qe_file_list)

task.run(args.sensor_id, qe_files, args.pd_ratio_file,
         args.mask_files(qe_files[0]), args.system_gains(),
         medians_file=args.medians_file, vendor_data=args.vendor_data,
         correction_image=args.correction_image)
