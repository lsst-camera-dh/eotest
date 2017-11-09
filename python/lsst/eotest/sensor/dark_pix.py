import sys
import traceback
import numpy as np
import argparse
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.eotest.image_utils as imutils


class DarkPix(object):
    def __init__(self, percent):
        self.percent = percent

    def __call__(self, fitsfile, amps):

        tot_dark_ccd = 0
        tot_dark_per_amp = []
        pix_per_amp = []
        col_per_amp = []

        for amp in amps:
            try:
                tot, pixels = self.dark_pix(fitsfile, imutils.dm_hdu(amp))
                tot_dark_ccd += tot
                tot_dark_per_amp.append(tot)
                pix_per_amp.append(pixels)
            except:
                print "Failed dark pixels for hdu ", amp, " ", \
                    imutils.dm_hdu(amp)
                traceback.print_exc(file=sys.stdout)
                continue

        return tot_dark_ccd, tot_dark_per_amp, pix_per_amp, col_per_amp

    def dark_pix(self, infile, hdu):
        """ List pixels with counts less than a specified percentage of the
            median flux, for the specified amps. """

        #read in and trim image area
        im = imutils.unbias_and_trim(afwImage.ImageF(infile, hdu))

        #find median of image
        median = afwMath.makeStatistics(im, afwMath.MEDIAN).getValue()
        thresh = median*self.percent/100.0

        #find pixels less than _ percent of median
        imarr = im.getArray()
        darkpix = np.where(imarr <= thresh)

        #turn x,y into a list of pixel coords
        pixlist = np.transpose(darkpix)

        return len(pixlist), pixlist


def run_dark_pix(fitsfile, percent=80, amps=None, verbose=False):
    """ Given an input FITS file, find bright pixels."""
    if amps is None:
        amps = imutils.allAmps(fitsfile)
    try:
        dp = DarkPix(percent)
        tot_dark_pixels, tot_per_amp, pix_per_amp, col_per_amp = dp(fitsfile, amps)
    except:
        traceback.print_exc(file=sys.stdout)

    if verbose:
        for ind, amp in enumerate(amps):
            print "Amp: ", amp, " ", tot_per_amp[ind], " Dark Pixels"
            print pix_per_amp[ind]
        print 'Total CCD Dark Pixels: ', tot_dark_pixels

    return tot_dark_pixels, tot_per_amp, pix_per_amp, col_per_amp


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Find the locations of dark pixels.')
    parser.add_argument('-i', '--infile', help="path to input image file",
                        type=str)
    parser.add_argument('-p', '--percent',
                        help="Percentage of median to use as threshold",
                        type=float, default=80.0)
    parser.add_argument('-a', '--amps',
                        help="amps to be analyzed, separated by a space",
                        type=int, nargs='+', default=range(1, 17))
    parser.add_argument('-v', '--verbose', help="turn verbosity on",
                        action='store_true', default=False)
    args = parser.parse_args()

    tot, tot_per_amp, pix_per_amp, col_per_amp = \
        run_dark_pix(args.infile, args.percent, args.amps, args.verbose)
