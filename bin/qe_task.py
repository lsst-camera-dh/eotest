#!/usr/bin/env python

"""
@brief Compute QE curves from the wavelength scan dataset.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import lsst.eotest.sensor as sensorTest
import lsst.eotest.sensor.qe as qe

if __name__ == '__main__':
    parser = sensorTest.TaskParser('Compute QE curves')
    parser.add_argument('-f', '--qe_files', type=str,
                        help='wavelength scan file pattern')
    parser.add_argument('-F', '--qe_file_list', type=str,
                        help='list of wavelength scan files')
    parser.add_argument('-q', '--qe_medians_file', default=None, type=str,
                        help='file of median image values from wl scan dataset')
    parser.add_argument('-c', '--ccd_cal_file', type=str,
                        help='calibration file for photodiode at CCD location')
    parser.add_argument('-i', '--int_sph_cal_file', type=str,
                        help='calibration file for photodiode at integrating sphere')
    parser.add_argument('-w', '--wavelength_scan_file', type=str,
                        help='uncorrected wavelength scan file')
    args = parser.parse_args()

    gains = args.system_gains()
    sensor = args.sensor()
    sensor_id = args.sensor_id
    
    if args.qe_medians_file is None:
        medians_file = os.path.join(args.output_dir, '%s_QE.txt' % sensor_id)
    else:
        medians_file = args.qe_medians_file

    fits_outfile = os.path.join(args.output_dir, '%s_QE.fits' % sensor_id)

    ccd_cal_file = args.ccd_cal_file
    sph_cal_file = args.int_sph_cal_file
    wlscan_file = args.wavelength_scan_file

    infiles = args.files(args.qe_files, args.qe_file_list)

    qe_data = qe.QE_Data()
    if args.qe_medians_file is None:
        qe_data.calculate_medians(infiles, medians_file,
                                  mask_files=args.mask_files(),
                                  clobber=True)
    qe_data.read_medians(medians_file)
    qe_data.calculate_QE(ccd_cal_file, sph_cal_file, wlscan_file, gains)
    qe_data.write_fits_tables(fits_outfile)
