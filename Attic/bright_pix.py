"""
@brief Find bright pixels and bright columns above a 
threshold = mean + nsig*sigma
"""
import argparse
import numpy as np
import lsst.afw.image as afwImage
import lsst.daf.base as dafBase
import image_utils
import sys, traceback
#from image_utils import unbias_and_trim
from simulation.sim_tools import SegmentExposure, writeFits

class BrightPix(object):
    def __init__(self, nsig=5):
        self.nsig = nsig

    def __call__(self, fitsfile, amps):
        """ Iterate over requested amps and find bright pixels and columns. """

        tot_bright_ccd = 0
        tot_bright_per_amp = []
        pix_per_amp = []
        col_per_amp = []

        for amp in amps:
            try:
                tot, pixels, cols = self.bright_pix(fitsfile, \
                                                    image_utils.dm_hdu(amp))
                tot_bright_ccd += tot
                tot_bright_per_amp.append(tot)
                pix_per_amp.append(pixels)
                col_per_amp.append(cols)
            except:
                print "Failed bright pixel for hdu ", amp, " ", \
                      image_utils.dm_hdu(amp)
                traceback.print_exc(file=sys.stdout)
                continue

        return tot_bright_ccd, tot_bright_per_amp, pix_per_amp, col_per_amp

    def bright_pix(self, fitsfile, hdu):
        """ Does the work of finding the bright pixels and columns. """

        im = image_utils.unbias_and_trim(afwImage.ImageF(fitsfile, hdu))
        imarr = im.getArray()
        mean = np.mean(imarr)
        sigma = np.std(imarr)
        threshold = self.nsig*sigma + mean

        # Find bright pixels.
        pixels = np.where(imarr > threshold)

        # Find bright columns.
        col_means = [np.mean(imarr[:, i]) for i in range(im.getWidth())]
        columns = np.where(col_means > threshold)

        # Weed out bright pixels that are already in bright columns or rows.
        indx = [i for i in range(len(pixels[1])) 
                if pixels[1][i] not in columns[0]]

        pixels = (pixels[0][indx], pixels[1][indx])
        tup = zip(pixels[1], pixels[0])
        sorted_tup = sorted(tup)
        return len(sorted_tup), sorted_tup, columns


def run_bright_pix(fitsfile, amps=image_utils.allAmps, verbose=False):
    """ Given an input FITS file, find bright pixels."""

    try:
        bp = BrightPix()
        tot_bright_pixels, tot_per_amp, pix_per_amp, col_per_amp = bp(fitsfile, amps)
    except:
        traceback.print_exc(file=sys.stdout)

    if verbose:
        for ind, amp in enumerate(amps):
            print "Amp: ", amp, " ", tot_per_amp[ind], " Bright Pixels"
            print pix_per_amp[ind]
        print 'Total CCD Bright Pixels: ', tot_bright_pixels

    return tot_bright_pixels, tot_per_amp, pix_per_amp, col_per_amp


def write_test_image(outfile, nhdus=16, verbose=True):
    if verbose:
        print "simulating bright pixels:", outfile
    segments = []
    for hdu in range(nhdus):
        seg = SegmentExposure()
        seg.add_bias(1e4, 10)
        seg.add_dark_current(300)
        seg.expose_flat(200)
        cols = seg.add_bright_cols(ncols=1, nsig=10)
        pix = seg.add_bright_pix(npix=100, nsig=10)
        segments.append(seg)
    writeFits(segments, outfile)

def run_test():
    """ Generate a test file and find bright pixels. """

    fitsfile = 'test_image.fits'
    write_test_image(fitsfile)
    tot_bright_pix, tot_per_amp, tup_per_amp, col_per_amp = \
        run_bright_pix(fitsfile, verbose=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Find bright pixels.')
    parser.add_argument('-i', '--infile', \
                        help="image file to be used for analysis", type=str)
    parser.add_argument('-a', '--amps', \
                        help="amps to be analyzed, separated by a space", \
                        type=int, nargs = '+', default=range(1, 17))
    parser.add_argument('-t', '--test', help="run test only", \
                        action='store_true', default=False)
    parser.add_argument('-v', '--verbose', help="turn verbosity on", \
                        action='store_true', default=False)
    args = parser.parse_args()

    if args.test:
        run_test()
    else:
        tot, tot_per_amp, pix_per_amp, col_per_amp = \
            run_bright_pix(args.infile, args.amps, args.verbose)


