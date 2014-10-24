#!/usr/bin/env python

import os
import glob
import numpy as np
import pyfits
import lsst.eotest.image_utils as imutils

def eotest_check_input_data(rootdir='.', use_baselined=True):
    """
    Check electro-optical datasets according to relevant specs in
    LCA-128-E, LCA-10103-A, and LCA-10140-A?. Not all specs are tested.

    rootdir is assumed to point to a directory that ends in the sensor
    ID, i.e., it is of the form "<...>/NNN-MM".  Checks are performed
    only for files and directories below rootdir.
    """
    errors = {'dirs' : [], 'temperature' : [], 'lambda' : []}
    full_path = lambda x : os.path.join(rootdir, x)
    if use_baselined:  # Use the baselined version in docushare
        print "Assuming baselined version of LCA-10140-A (from docushare)"
        test_types = ('dark', 'fe55', 'flat', 'lambda', 'spot',
                      'superflat_500', 'trap')
        file_pattern = os.path.join('*', '*.fits')
    else:           # Use the version in the most recent revision (2014-04-11)
        print "Assuming 2014-04-11 revision of LCA-10140-A"
        test_types = ('dark', 'fe55', 'flat', 'lambda', 'spot', 
                      'sflat_500', 'trap')
        file_pattern = os.path.join('*', '*', '*.fits')
    files = {}
    #
    # Check basic directory structure and glob for filenames as
    # specified in LCA-10140-A.
    #
    for test in test_types:
        target = full_path(test)
        if not os.path.isdir(target):
            what = "Expected test type subdir %(target)s not found" % locals()
            errors['dirs'].append(what)
        glob_target = os.path.join(rootdir, test, file_pattern)
        files[test] = glob.glob(glob_target)
        if len(files[test]) == 0:
            what = "No files found for test type %(test)s in %(glob_target)s" \
                % locals()
            errors['dirs'].append(what)
    #
    # Check temperature ranges of datasets need for read noise (fe55),
    # crosstalk (spot), dark current (dark), and QE (lambda) as
    # specified in LCA-128-E.
    #
    for test in ('fe55', 'spot', 'dark', 'lambda'):
        for infile in files[test]:
            try:
                imutils.check_temperatures([infile], 0.1, setpoint=-95)
            except RuntimeError, eObj:
                errors['temperature'].append(eObj.message)
    #
    # Check for required wavelengths for QE and PRNU as specified in
    # LCA-128-E.
    #
    required_wls = (330, 350, 370, 450, 500, 620, 750, 870, 1000)
    acquired_wls = []
    for item in files['lambda']:
        wl = int(np.round(pyfits.open(item)[0].header['MONOWL']))
        acquired_wls.append(wl)
    for target_wl in required_wls:
        if target_wl not in acquired_wls:
            what = 'Flat at %(target_wl)s nm missing from lambda dataset' \
                % locals()
            errors['lambda'].append(what)
    #
    # Print summary of failures.
    #
    for test in errors:
        if errors[test]:
            print "%i error(s) of type %s:" % (len(errors[test]), test)
            for message in errors[test]:
                print "  ", message
            print

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Check EO Test dataset compliance')
    parser.add_argument('sensor_dir', help='Root directory for data from a specific sensor.  Should point to a directory of the form <...>/NNN-MM')
    parser.add_argument('-b', '--use_baselined', action='store_true', 
                        help='Flag to use baselined version of LCA-10140 (in docushare)')
    args = parser.parse_args()
    eotest_check_input_data(args.sensor_dir, args.use_baselined)
