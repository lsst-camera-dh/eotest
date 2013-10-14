#!/usr/bin/env python

"""
@brief Pixel response non-uniformity.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import lsst.test_scripts.sensor as sensorTest
#from MaskedCCD import MaskedCCD, Metadata
#from pipeline.TaskParser import TaskParser
#from prnu import prnu

if __name__ == '__main__':
    parser = sensorTest.TaskParser('Compute pixel response non-uniformity')
    parser.add_argument('-f', '--file_pattern', type=str,
                        help='file pattern for PRNU input files')
    parser.add_argument('-F', '--prnu_file_list', type=str,
                        help='list of files to be used in PRNU calculations')
    parser.add_argument('-c', '--correction_image', type=str,
                        help='correction image for illumination non-uniformity',
                        default=None)
    args = parser.parse_args()
    sensor = args.sensor()
    sensor_id = args.sensor_id
    mask_files = args.mask_files()
    gains = args.system_gains()
    correction_image = args.correction_image
    outfile = os.path.join(args.output_dir, '%s_prnu_values.txt' % sensor_id)

    infiles = args.files(args.file_pattern, args.prnu_file_list)

    output = open(outfile, 'w')
    output.write("wavelength (nm)   pixel_stdev   pixel_median\n")
    for infile in infiles:
        md = sensorTest.Metadata(infile, 1)
        wl = md.get('MONOWL')
        if int(wl) in (350, 450, 500, 620, 750, 870, 1000):
            print "  working on", wl
            pix_stdev, pix_median = sensorTest.prnu(infile, mask_files, gains,
                                                    correction_image=correction_image)
            output.write("%6.1f  %12.4e  %12.4e\n"
                         % (wl, pix_stdev, pix_median))
            output.flush()
    output.close()
