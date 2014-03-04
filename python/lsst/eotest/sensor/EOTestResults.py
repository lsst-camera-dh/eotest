import os
import numpy as np
import pyfits

_namps = 16

class EOTestResults(object):
    def __init__(self, infile):
        self.outfile = infile
        self.extname = 'AMPLIFIER_RESULTS'
        if not os.path.isfile(infile):
            self._createFitsObject()
        else:
            self.output = pyfits.open(infile)
    def _createFitsObject(self):
        self.output = pyfits.HDUList()
        self.output.append(pyfits.PrimaryHDU())
        self.colnames = ["AMP", "GAIN", "READ_NOISE", "FULL_WELL",
                         "CTI_SERIAL", "CTI_PARALLEL", "DARK_CURRENT_95",
                         "NUM_BRIGHT_PIXELS"]
        formats = "IEEEEEEI"
        my_types = dict((("I", np.int), ("E", np.float)))
        columns = [np.zeros(_namps, dtype=my_types[fmt]) for fmt in formats]
        units = ["None", "Ne/DN", "rms e-/pixel", "e-/pixel",
                 "None", "None", "e-/s/pixel", "None"]
        fits_cols = [pyfits.Column(name=self.colnames[i], format=formats[i],
                                   unit=units[i], array=columns[i])
                     for i in range(len(self.colnames))]
        self.output.append(pyfits.new_table(fits_cols))
        self.output[-1].name = self.extname
        for amp in range(1, _namps+1):
            self.add_seg_result(amp, 'AMP', amp)
    def __getitem__(self, column):
        return self.output[self.extname].data.field(column)
    def add_seg_result(self, amp, column, value):
        self.output[self.extname].data.field(column)[amp-1] = value
    def write(self, outfile=None, clobber=True):
        if outfile is None:
            outfile = self.outfile
        self.output.writeto(outfile, clobber=clobber)

if __name__ == '__main__':
    outfile = 'foo.fits'
    foo = EOTestResults(outfile)
    for amp in range(1, 17):
        foo.add_seg_result(amp, 'GAIN', 5)
    foo.write()

    bar = EOTestResults(outfile)
    print bar['GAIN']

