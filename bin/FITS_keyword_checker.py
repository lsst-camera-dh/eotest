#!/usr/bin/env python
"""
@brief Check that FITS files have all of the header keywords specified in
LCA-10140.
"""

import sys
import argparse
import lsst.eotest.sensor as sensorTest

parser = argparse.ArgumentParser(description='Check FITS header keywords for sensor data files.')

parser.add_argument('infile', help='CCD image file to be checked')
parser.add_argument('-u', '--used_keywords', action='store_true',
                    default=False,
                    help='check only keywords used by test scripts')
parser.add_argument('-v', '--verbose', action='store_true', default=False,
                    help='verbosity (print results to screen)')

args = parser.parse_args()

if args.used_keywords:
    template = sensorTest.fits_headers.template_used_file
else:
    template = sensorTest.fits_headers.template_file
    
missing = sensorTest.fits_headers.check_keywords(args.infile,
                                                 template=template,
                                                 verbose=args.verbose)

# Use the total number of missing keywords as the return code.
retcode = sum(len(x) for x in missing.values())
sys.exit(retcode)
