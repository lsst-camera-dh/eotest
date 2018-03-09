"""
@brief Class to contain results from EO testing and save in a FITS file
as a binary table.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
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
    def _createFitsObject(self):
        self.output = fits.HDUList()
        self.output.append(fits.PrimaryHDU())
        self.colnames = ["AMP", "GAIN", "GAIN_ERROR", "READ_NOISE", "FULL_WELL",
                         "CTI_HIGH_SERIAL", "CTI_HIGH_PARALLEL",
                         "CTI_LOW_SERIAL", "CTI_LOW_PARALLEL",
                         "DARK_CURRENT_95", "NUM_BRIGHT_PIXELS", "NUM_TRAPS"]
        formats = "IEEEEEEEEEII"
        my_types = dict((("I", np.int), ("E", np.float)))
        columns = [np.zeros(self.namps, dtype=my_types[fmt]) for fmt in formats]
        units = ["None", "Ne/DN", "Ne/DN", "rms e-/pixel", "e-/pixel",
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
        _types = dict(((int, 'I'), (float, 'E'), (np.float64, 'E')))
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
        fitsWriteto(self.output, outfile, clobber=clobber)
    def defects_fractions(self, col_len=None, total_pixels=None):
        """
        Sum the bad pixel contributions from the various defect types
        and express as a fraction per segment.  This is based on POC's
        algorithm as given in LSSTTD-1255.

        Returns
        -------
        np.array: array of defects fractions for all of the amplifiers
           in the sensor.
        """
        vendor_col_len = {'ITL': 2000, 'E2V': 2002}
        vendor_total_pixels = {'ITL': 1024000, 'E2V': 1025024}
        # Extract the vendor from the filename.  TODO: extract this
        # from FITS keywords instead.
        vendor = self.infile[:3]
        if col_len is None:
            col_len = vendor_col_len[vendor]
        if total_pixels is None:
            total_pixels = vendor_total_pixels[vendor]
        bad_pix = (self['NUM_BRIGHT_PIXELS'][amp-1]
                   + self['NUM_DARK_PIXELS'][amp-1]
                   + col_len*(self['NUM_BRIGHT_COLUMNS']
                              + self['NUM_DARK_COLUMNS']))
        return float(bad_pix)/total_pixels
    def sensor_grade(self):
        pass


if __name__ == '__main__':
    outfile = 'foo.fits'
    foo = EOTestResults(outfile)
    print foo.colnames
    for amp in range(1, 17):
        foo.add_seg_result(amp, 'GAIN', 5)
        foo.add_seg_result(amp, 'NEW_INT_COLUMN', 2)
    foo.append_column('NEW_FLOAT_COLUMN', type(3.))
    foo.write()

    bar = EOTestResults(outfile)
    print bar['GAIN']
    print bar['NEW_INT_COLUMN']
    print bar['NEW_FLOAT_COLUMN']
