"""
@brief Class to compute medians per amplifier for a set of exposures taken
as a function of wavelength.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import glob
import numpy as np
import lsst.afw.math as afwMath
import image_utils as imutils
from MaskedCCD import MaskedCCD, Metadata

class QE_Data(object):
    def __init__(self, results_file, pattern=None, verbose=False):
        self.results_file = results_file
        self.verbose = verbose
        if pattern is not None:
            self._calculate_medians(pattern)
        else:
            self._read_medians()
    def _calculate_medians(self, pattern):
        files = glob.glob(pattern)
        files.sort()
        self.medians = dict([(amp, []) for amp in imutils.allAmps])
        self.wl = []

        output = open(self.results_file, 'w')
        for item in files[:10]:
            if self.verbose:
                print "processing", item
            ccd = MaskedCCD(item)
            md = Metadata(item, 1)
            wl.append(md.get('MONO.WAVELENG'))
            output.write('%.3f' % wl[-1])
            for amp in imutils.allAmps:
                im = imutils.unbias_and_trim(ccd[amp])
                value = afwMath.makeStatistics(im, afwMath.MEDIAN).getValue()
                medians[amp].append(value)
                output.write('  %.5f' % medians[amp][-1])
            output.write('\n')
            output.flush()
        output.close()
        self.wl = np.array(self.wl)
    def _read_medians(self):
        data = np.recfromtxt(self.results_file)
        data = data.transpose()
        self.wl = data[0]
        self.medians = dict([(amp, col) for amp, col
                             in zip(imutils.allAmps, data[1:])])

if __name__ == '__main__':
    outfile = 'med_vs_wl.txt'
    pattern = '/nfs/farm/g/lsst/u1/testData/HarvardData/112-01/final/bss50/qe'
    
    qe_data = QE_Data(outfile)
