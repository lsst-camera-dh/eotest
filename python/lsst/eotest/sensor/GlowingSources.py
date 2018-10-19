"""
@brief Find extended bright sources in a dark image.
"""
import os
import numpy as np
import astropy.io.fits as fits
from lsst.eotest.fitsTools import fitsWriteto, fitsTableFactory

import lsst.afw.detection as afwDetect
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath

import lsst.eotest.image_utils as imutils
#from .MaskedCCD import MaskedCCD
#from .AmplifierGeometry import makeAmplifierGeometry
#from .fits_headers import fits_headers
from MaskedCCD import MaskedCCD
from AmplifierGeometry import makeAmplifierGeometry
from fits_headers import fits_headers

def prect_values(peak, imarr, ixm=10, ixp=10, iym=10, iyp=10):
    """Make a postage stamp of a detected source with given dimensions.
    """
    xpeak, ypeak = peak.getIx(), peak.getIy()
    # nb. imarr includes overscan
    yimsiz,ximsiz = imarr.shape
    # if we are too close to an edge, just return all 0's
    if ypeak-iym<0 or xpeak-ixm<0 or ypeak+iyp+1>yimsiz or xpeak+ixp+1>ximsiz:
        prect_data = np.zeros([iyp+iym+1, ixp+ixm+1])
    else:
        # store a region around peak pixel [-3,21] in x and [-3,3] in y
        prect_data = imarr[ypeak-iym:ypeak+iyp+1, xpeak-ixm:xpeak+ixp+1]

    # data is stored as a normal 2-d array and then flattened
    prect_data_flat = prect_data.flatten()
    if prect_data_flat.shape[0] == 0:
        pdb.set_trace()
    return prect_data_flat

class GlowingSources(object):
    """
    Find extended bright sources in a dark image for signals greater than
    nsig above the median and extended over at least npix_min.

    Detected sources are saved in a FITs table.
    """
    def __init__(self, nsig=10, npix_min=20, outfile=None):

        self.nsig = nsig
        self.npix_min = npix_min

        self.dn_fp = []
        self.amp = []
        self.amp_set = set()

        ## Add an output file
        self.outfile = outfile
        if outfile is None:
            self.output = fits.HDUList()
            self.output.append(fits.PrimaryHDU())
        else:
            self.output = fits.open(self.outfile)
     
    def find(self, ccd, amp, oscan_fit_order=1):
        """
        Find and return extended glowing sources.
        """
        image = ccd.bias_subtracted_image(amp, fit_order=oscan_fit_order)
        imarr = image.getImage().getArray()

        flags = afwMath.MEDIAN | afwMath.STDEVCLIP
        statistics = afwMath.makeStatistics(image, flags, ccd.stat_ctrl)
        median = statistics.getValue(afwMath.MEDIAN)
        stdev = statistics.getValue(afwMath.STDEVCLIP)

        threshold = afwDetect.Threshold(median + self.nsig*stdev)
        fpset = afwDetect.FootprintSet(image.getImage(), threshold, 
                                       self.npix_min)

        x0, y0 = [], []
        dn_fp = []
        prect_data = []
        num_fp = 0
        for fp in fpset.getFootprints():
            num_fp += 1
            spans = fp.getSpans()
            positions = []
            zvals = []
            peak = [pk for pk in fp.getPeaks()][0]
            dn_sum = 0
            for span in spans:
                y = span.getY()
                for x in range(span.getX0(), span.getX1() + 1):
                    positions.append((x, y))
                    zvals.append(imarr[y][x])
                    dn_sum += imarr[y][x]
            x0.append(peak.getIx())
            y0.append(peak.getIy())
            dn_fp.append(dn_sum)
            prect_data_row = prect_values(peak,imarr)
            prect_data.append(prect_data_row)

        self._save_ext_data(amp, x0, y0, dn_fp, np.array(prect_data))
        self.amp_set.add(amp)
        self.dn_fp.extend(dn_fp)
        self.amp.extend(np.ones(len(x0))*amp)

    def _save_ext_data(self, amp, x0, y0, dn_fp, prect_data):
        """Save an FITs table extension for each amplifier."""

        print amp, x0, y0, dn_fp

        extname = 'Amp%02i' % amp
        try:
            table_hdu = self.output[extname]
            row0 = table_hdu.header['NAXIS2']
            nrows = row0 + len(x0)
            table_hdu = fitsTableFactory(table_hdu.data, nrows=nrows)
            for i in range(len(x0)):
                row = i + row0
                table_hdu.data[row]['AMPLIFIER'] = amp
                table_hdu.data[row]['XPOS'] = x0[i]
                table_hdu.data[row]['YPOS'] = y0[i]
                table_hdu.data[row]['DN'] = dn[i]
                table_hdu.data[row]['PRECT_DATA'] = prect_data[i]
            table_hdu.name = extname
            self.output[extname] = table_hdu
        except KeyError:

            colnames = ['AMPLIFIER', 'XPOS', 'YPOS', 'DN_FP', 'PRECT_DATA']
            columns = [np.ones(len(x0))*amp, np.array(x0), np.array(y0),
                       np.array(dn_fp), np.array(prect_data)]
            formats = ['I'] + ['E']*3 + ['441E']
            units = ['None', 'pixel', 'pixel', 'ADU', 'ADU']
            fits_cols = lambda coldata: [fits.Column(name=colname, 
                                                     format=format,
                                                     unit=unit,
                                                     array=column)
                                         for colname, format, unit, column
                                         in coldata]
            self.output.append(fitsTableFactory(fits_cols(zip(colnames,
                                                             formats,
                                                             units,
                                                             columns))))
            self.output[-1].name = extname

    def write_results(self, outfile='glowing_sources_params.fits'):
        """Write results to the corresponding output FITs file."""
        self.output[0].header['NAMPS'] = len(self.amp_set)
        fitsWriteto(self.output, outfile, overwrite=True, checksum=True)

if __name__ == '__main__':
    import argparse
 
    parser = argparse.ArgumentParser()
    parser.add_argument('dark_file', 
                        help='Dark file to use for bad pixel identification.')
    parser.add_argument('--nsig', '-n', type=int, default=10)
    parser.add_argument('--min_pixels', '-p', type=int, default=20)
    args = parser.parse_args()

    dark_file = args.dark_file
    nsig = args.nsig
    min_pixels = args.min_pixels
    gain = 0.7

    ccd = MaskedCCD(dark_file)
    glows = GlowingSources(nsig=nsig, npix_min=min_pixels)
    for amp in ccd:
        glows.find(ccd, amp)
    glows.write_results()
