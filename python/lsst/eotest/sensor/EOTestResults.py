"""
@brief Class to contain results from EO testing and save in a FITS file
as a binary table.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
from __future__ import print_function
import os
from collections import OrderedDict
import numpy as np
import astropy.io.fits as fits
from lsst.eotest.fitsTools import fitsTableFactory, fitsWriteto


class EOTestResults(object):
    """
    This class saves EO test results per segment/amplifier in a FITS
    binary table.  The data to be collected are specified in LCA-10301-A.
    """

    def __init__(self, infile, namps=16):
        self.infile = infile
        self.namps = namps
        self.extname = 'AMPLIFIER_RESULTS'
        if not os.path.isfile(infile):
            self._createFitsObject()
        else:
            self.output = fits.open(infile)
            self.colnames = self.output[self.extname].data.names

    def __del__(self):
        self.output.close()

    def _createFitsObject(self):
        self.output = fits.HDUList()
        self.output.append(fits.PrimaryHDU())
        self.colnames = ["AMP", "GAIN", "GAIN_ERROR", "READ_NOISE", "FULL_WELL",
                         "MAX_FRAC_DEV",
                         "CTI_HIGH_SERIAL", "CTI_HIGH_PARALLEL",
                         "CTI_LOW_SERIAL", "CTI_LOW_PARALLEL",
                         "DARK_CURRENT_95", "NUM_BRIGHT_PIXELS", "NUM_TRAPS"]
        formats = "IEEEEEEEEEEJJ"
        my_types = dict((("I", np.int), ("J", np.int), ("E", np.float)))
        columns = [np.zeros(self.namps, dtype=my_types[fmt]) for fmt in formats]
        units = ["None", "Ne/DN", "Ne/DN", "rms e-/pixel", "e-/pixel", "None",
                 "None", "None", "None", "None", "e-/s/pixel", "None", "None"]
        fits_cols = [fits.Column(name=self.colnames[i], format=formats[i],
                                 unit=units[i], array=columns[i])
                     for i in range(len(self.colnames))]
        self.output.append(fitsTableFactory(fits_cols))
        self.output[-1].name = self.extname
        for amp in range(1, self.namps+1):
            self.add_seg_result(amp, 'AMP', amp)

    def __getitem__(self, column):
        try:
            return self.output[self.extname].data.field(column)
        except:
            return self.output[column]

    def append_column(self, colname, dtype=np.float, unit='None', column=None):
        """
        Append a new column of amplifier data to the AMPLIFIER_RESULTS table.
        """
        if colname in self.colnames:
            return
        _types = dict(((int, 'J'), (float, 'E'), (np.float32, 'E'),
                       (np.float64, 'E')))
        if column is None:
            column = np.zeros(self.namps, dtype=dtype)
        new_cols = fits.ColDefs([fits.Column(name=colname,
                                             format=_types[dtype],
                                             unit=unit, array=column)])
        new_hdu = fitsTableFactory(self.output[self.extname].data.columns
                                   + new_cols)
        new_hdu.name = self.extname
        self.output[self.extname] = new_hdu
        self.colnames.append(colname)

    def add_seg_result(self, amp, column, value):
        """
        Add the results for a given amplifier segment and column.
        """
        if column not in self.colnames:
            self.append_column(column, type(value))
        self.output[self.extname].data.field(column)[amp-1] = value

    def add_ccd_result(self, keyword, value):
        """
        Add CCD-wide key/value pair to the primary HDU of the output.
        """
        self.output[0].header[keyword] = value

    def write(self, outfile=None, clobber=True):
        """
        Write or update the output file.
        """
        if outfile is None:
            outfile = self.infile
        fitsWriteto(self.output, outfile, overwrite=clobber)

    def defect_fractions(self, col_len=None, total_pixels=None):
        """
        Sum the bad pixel contributions from the various defect types
        and express as a fraction per segment.  This is based on POC's
        algorithm as given in LSSTTD-1255.

        Returns
        -------
        np.array: array of defect fractions for all of the amplifiers
           in the sensor.
        """
        vendor_col_len = {'ITL': 2000, 'E2V': 2002}
        vendor_total_pixels = {'ITL': 1024000, 'E2V': 1025024}
        # Get vendor from FITS header, if set; otherwise, extract the
        # from the filename, assuming the standard naming convention.
        try:
            vendor = self.output[0].header['CCD_MANU']
        except KeyError:
            vendor = self.infile[:3]
        if col_len is None:
            col_len = vendor_col_len[vendor]
        if total_pixels is None:
            total_pixels = vendor_total_pixels[vendor]
        bad_pix = (self['NUM_BRIGHT_PIXELS']
                   + self['NUM_DARK_PIXELS']
                   + col_len*(self['NUM_BRIGHT_COLUMNS']
                              + self['NUM_DARK_COLUMNS']))
        return bad_pix/float(total_pixels)

    @staticmethod
    def bias_offset_grade(bnl_bias_stats_file=None):
        if bnl_bias_stats_file is not None:
            with open(bnl_bias_stats_file) as input_:
                bias_levels = eval(input_.readline().strip().split(':')[1])
                med1_8 = np.median(bias_levels[:8])
                med9_16 = np.median(bias_levels[8:])
                if (max(med1_8, med9_16) + 3000.)/2. > 6000.:
                    return '-HI\_BIAS'
        return ''

    def sensor_stats(self, bnl_bias_stats_file=None):
        """
        Return the sensor statistics based on the EO test results.  This
        applies POC's sensor grade criteria as given in LSSTTD-1255.

        Returns
        -------
        str: "SCIENCE", "RESERVE", or "ENGIN."
        """
        rn = self['READ_NOISE']
        ppm = 1e6
        cti_ls = ppm*self['CTI_LOW_SERIAL']
        cti_hs = ppm*self['CTI_HIGH_SERIAL']
        cti_lp = ppm*self['CTI_LOW_PARALLEL']
        cti_hp = ppm*self['CTI_HIGH_PARALLEL']
        defects = self.defect_fractions().mean()
        num_bright_cols = self['NUM_BRIGHT_COLUMNS'].sum()
        if (rn.max() < 10 and
            (rn > 8.5).sum() < 3 and
            cti_ls.max() < 9 and
            cti_lp.max() < 7 and
            cti_hs.max() < 9 and
            cti_hp.max() < 5 and
            defects < 0.0049 and
            num_bright_cols < 5):
            GRADE = "SCIENCE" + self.bias_offset_grade(bnl_bias_stats_file)
        elif (rn.max() < 14 and
              (rn > 11).sum() < 4 and
              cti_ls.max() < 18 and
              (cti_ls > 10).sum() < 4 and
              cti_lp.max() < 10 and
              cti_hp.max() < 10 and
              (cti_lp > 10).sum() < 4):
            GRADE = "RESERVE"
        else:
            GRADE = "ENGIN."

        stats = OrderedDict()
        stats['GRADE'] = GRADE
        stats['max rn'] = rn.max()
        stats['\# rn>8'] = (rn > 8).sum()
        stats['max SCTI'] = max(cti_ls.max(), cti_hs.max())
        stats['\# SCTI>5'] = (cti_ls > 5).sum()
        stats['max PCTI'] = max(cti_lp.max(), cti_hp.max())
        stats['\# PCTI>3'] = (cti_lp > 3).sum()
        stats['\% defects'] = 100*defects
        stats['\# bright cols'] = num_bright_cols
        return stats


if __name__ == '__main__':
    outfile = 'foo.fits'
    foo = EOTestResults(outfile)
    print(foo.colnames)
    for amp in range(1, 17):
        foo.add_seg_result(amp, 'GAIN', 5)
        foo.add_seg_result(amp, 'NEW_INT_COLUMN', 2)
    foo.append_column('NEW_FLOAT_COLUMN', type(3.))
    foo.write()

    bar = EOTestResults(outfile)
    print(bar['GAIN'])
    print(bar['NEW_INT_COLUMN'])
    print(bar['NEW_FLOAT_COLUMN'])
