"""
@brief Find traps from a pocket-pumped CCD frame.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
from __future__ import absolute_import
import os
import numpy as np
import astropy.io.fits as fits
from lsst.eotest.fitsTools import fitsTableFactory, fitsWriteto
import lsst.afw.detection as afwDetect
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.eotest.image_utils as imutils
from .MaskedCCD import MaskedCCD
from .TrapFinder import TrapFinder


class Traps(dict):
    """
    Find traps from a pocket-pumped image, listing the trap x, y
    positions and size by amplifier.  Provide a method to write the
    results to a FITS file as a binary table.
    """

    def __init__(self, ccd, gains, cycles=100, C2_thresh=10,
                 C3_thresh=1, nx=10, ny=10,
                 edge_rolloff=10, amps=None):
        super(Traps, self).__init__()
        if amps is None:
            amps = list(ccd.keys())
        for amp in amps:
            self[amp] = []
            finder = TrapFinder(ccd, amp, C2_thresh=C2_thresh,
                                C3_thresh=C3_thresh, nx=nx, ny=ny,
                                edge_rolloff=edge_rolloff)
            results = finder.find()
            for ix, iy, c2, c3, a0, a1 in zip(*results):
                trap_size = max(a0, a1)*gains[amp]/cycles
                self[amp].append((ix, iy, trap_size, a0, a1))

    def write(self, outfile, overwrite=True):
        """
        Write the results as a FITS binary table.
        """
        nrows = sum([len(self[amp]) for amp in self])
        output = fits.HDUList()
        output.append(fits.PrimaryHDU())
        colnames = ['AMPLIFIER', 'XPOS', 'YPOS', 'TRAP_SIZE', 'A0', 'A1']
        formats = 'IIIIEE'
        units = ['None', 'pixel', 'pixel', 'electrons', 'ADU', 'ADU']
        columns = ([np.zeros(nrows, dtype=int)]*4
                   + [np.zeros(nrows, dtype=float)]*2)
        hdu = fitsTableFactory([fits.Column(name=colname, format=format,
                                            unit=unit, array=column)
                                for colname, format, unit, column in
                                zip(colnames, formats, units, columns)])
        hdu.name = 'TRAPS'
        output.append(hdu)
        row = 0
        for amp in self:
            for xpos, ypos, trap_size, a0, a1 in self[amp]:
                output['TRAPS'].data[row]['AMPLIFIER'] = amp
                output['TRAPS'].data[row]['XPOS'] = xpos
                output['TRAPS'].data[row]['YPOS'] = ypos
                output['TRAPS'].data[row]['TRAP_SIZE'] = trap_size
                output['TRAPS'].data[row]['A0'] = a0
                output['TRAPS'].data[row]['A1'] = a1
                row += 1
        fitsWriteto(output, outfile, overwrite=overwrite)


if __name__ == '__main__':
    gains = dict([(amp, 5) for amp in imutils.allAmps()])
    outfile = 'trap_list.fits'

#    infile = 'work/sensorData/000-00/trap/debug/000-00_trap_ppump_debug.fits'
    infile = '000-00_trap_ppump_debug.fits'
    mask_files = ()
    ccd = MaskedCCD(infile, mask_files=mask_files)

    my_traps = Traps(ccd, gains)
    my_traps.write(outfile)
