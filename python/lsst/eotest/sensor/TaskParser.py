"""
@brief Class to parse common command line options for pipeline tasks.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
from __future__ import print_function
from __future__ import absolute_import
from builtins import object
import os
import glob
import argparse
import astropy.io.fits as fits
import lsst.afw.image as afwImage
import lsst.eotest.image_utils as imutils
from .rolloff_mask import rolloff_mask
from .EOTestResults import EOTestResults
from lsst.eotest.database.SensorDb import SensorDb, NullDbObject
from lsst.eotest.database.SensorGains import SensorGains


class TaskNamespace(object):
    """
    Decorator class for argparse.Namespace.  This class provides
    functions for information derived from the Namespace attributes, such
    as lists of input files based on a file pattern or file_list file,
    a SensorDb object, system gains, etc..
    """

    def __init__(self, args):
        self.args = args

    def files(self, file_pattern, file_list):
        if file_list is not None:
            my_files = [x.strip() for x in open(file_list)
                        if x[0] != '#']
        elif file_pattern is not None:
            my_pattern = file_pattern.replace('\\', '')
            my_files = glob.glob(my_pattern)
            my_files.sort()
        else:
            raise RuntimeError('You must specify a file pattern or file list.')
        return my_files

    def bias_frame(self, outfile, overwrite=True):
        if (self.args.bias_frame_pattern is None and
                self.args.bias_frame_list is None):
            return None
        else:
            bias_files = self.files(self.args.bias_frame_pattern,
                                    self.args.bias_frame_list)
            if not bias_files:  # empty glob or file list
                return None
            imutils.fits_median_file(bias_files, outfile, clobber=True)
            return outfile

    def sensor(self):
        if self.args.db_credentials is not None:
            sensorDb = SensorDb(self.args.db_credentials)
            return sensorDb.getSensor(self.args.Vendor, self.args.sensor_id,
                                      add=True)
        else:
            return NullDbObject()

    def system_gains(self):
        if self.args.db_credentials is not None:
            return SensorGains(vendor=self.args.Vendor,
                               vendorId=self.args.sensor_id,
                               db_credentials=self.args.db_credentials)
        else:
            results = EOTestResults(self.args.gains)
            gains = results['GAIN']
            return dict([(amp, gains[amp-1]) for amp in imutils.allAmps()])

    def mask_files(self, infile=None):
        """
        Tuple of mask files to be used.  If infile is given and a
        rolloff defects mask is not found, then generate a roll-off
        mask file using infile to set the detector and amplifier
        geometry.
        """
        if self.args.verbose:
            print("mask file pattern:", self.args.mask_file_pattern)
        if self.args.mask_file_pattern is not None:
            my_mask_files = glob.glob(self.args.mask_file_pattern)
        else:
            my_mask_files = []
        # Loop through files and determine if edge rolloff mask is
        # included.
        have_rolloff_mask = False
        for mask_file in my_mask_files:
            hdr = fits.open(mask_file)[0].header
            if ('MASKTYPE' in list(hdr.keys()) and
                    hdr['MASKTYPE'] == 'ROLLOFF_DEFECTS'):
                have_rolloff_mask = True
        if infile is not None and not have_rolloff_mask:
            my_mask_files.append('edge_rollover_defect_mask.fits')
            if not os.path.isfile(my_mask_files[-1]):
                rolloff_mask(infile, my_mask_files[-1])
        if self.args.verbose:
            print("Using mask files:")
            for item in my_mask_files:
                print("  ", item)
        return my_mask_files

    def __getattr__(self, attrname):
        return getattr(self.args, attrname)


class TaskParser(argparse.ArgumentParser):
    """
    Subclass of argparse.ArgumentParser that adds default command line
    options for all pipeline tasks.
    """

    def __init__(self, description):
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
        argparse.ArgumentParser.__init__(self, description=description,
                                         formatter_class=formatter_class)
        self.add_argument('-d', '--db_credentials', type=str,
                          help='file containing database credentials')
        self.add_argument('-s', '--sensor_id', type=str,
                          help="sensor ID")
        self.add_argument('-b', '--bias_frame_pattern', type=str, default=None,
                          help='file pattern for bias frames')
        self.add_argument('-B', '--bias_frame_list', type=str, default=None,
                          help='file name of list of bias frames')
        self.add_argument('-V', '--Vendor', type=str,
                          help='CCD vendor (e.g., e2v, ITL)')
        self.add_argument('-m', '--mask_file_pattern', default=None, type=str,
                          help='file pattern for mask files')
        self.add_argument('-o', '--output_dir', type=str, default='.',
                          help="output directory")
        self.add_argument('-g', '--gains', type=str,
                          help='file of system gains')
        self.add_argument('-v', '--verbose', action='store_true', default=False,
                          help='turn verbosity on')
        self.add_argument('-r', '--results_file', type=str, default=None,
                          help='Results file for EO test parameters. Computed value if left at default of None: <SENSOR_ID>_eotest_results.fits')
        self.add_argument('--temp_set_point', type=float, default=-95.,
                          help='Expected temperture set point (in degrees C) for measurement')

    def parse_args(self):
        args = argparse.ArgumentParser.parse_args(self)
        if args.output_dir is not None:
            try:
                os.makedirs(args.output_dir)
            except OSError:
                pass
        return TaskNamespace(args)


if __name__ == '__main__':
    parser = TaskParser('test task')
    args = parser.parse_args()
