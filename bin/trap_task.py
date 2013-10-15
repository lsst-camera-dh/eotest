#!/usr/bin/env python

"""
@brief Task to find traps from pocket-pumped exposure.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import lsst.eotest.sensor as sensorTest

if __name__ == '__main__':
    parser = sensorTest.TaskParser('Find Charge Traps')
    parser.add_argument('-f', '--pocket_pumped_file', type=str,
                        help='Pocket pumped file')
    parser.add_argument('-O', '--output_file', type=str,
                        help='Output text file to contain trap locations and magnitudes')
    args = parser.parse_args()

    ccd = sensorTest.MaskedCCD(args.pocket_pumped_file,
                               mask_files=args.mask_files())

    if args.output_file is not None:
        outfile = os.path.join(args.output_dir, args.output_file)
    else:
        outfile = os.path.join(args.output_dir,
                               '%s_traps.txt' % args.sensor_id)
    sensorTest.traps(ccd, args.system_gains(), outfile=outfile)
