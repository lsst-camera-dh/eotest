"""
@brief Find bright pixels and bright columns above a 
threshold = mean + nsig*sigma
"""
import argparse
import numpy as np
import lsst.afw.image as afwImage

import image_utils
#from image_utils import unbias_and_trim
from simulation.sim_tools import SegmentExposure, writeFits

def bright_pix(infile, hdu=2, nsig=5):
    """ Does the work of finding the bright pixels and columns. """

    im = image_utils.unbias_and_trim(afwImage.ImageF(infile, hdu))
    imarr = im.getArray()

    mean = np.mean(imarr)
    sigma = np.std(imarr)
    threshold = nsig*sigma + mean

    # Find bright pixels.
    pixels = np.where(imarr > threshold)

    # Find bright columns.
    col_means = [np.mean(imarr[:, i]) for i in range(im.getWidth())]
    columns = np.where(col_means > threshold)

    # Weed out bright pixels that are already in bright columns or rows.
    indx = [i for i in range(len(pixels[1])) if pixels[1][i] not in columns]

    pixels = (pixels[0][indx], pixels[1][indx])

    return pixels, columns, im

def segment_results(amp, pixels, outfile=None, verbose=True):
    """ Print results for this segment/amplifier.

        A total count of bright pixels and table of locations
    """

    hdrStr = "Amp: " + str(amp) + " " + image_utils.hdu_dict[amp] \
              + " Bright Pixel Count: " + str(len(pixels[0]))
    if verbose:
        print hdrStr

    if outfile:
        outfile.write(hdrStr+"\n")

    tup = zip(pixels[0],pixels[1])
    sorted_tup = sorted(tup)
    for pair in sorted_tup:
        if verbose:
            print pair

        if outfile:
            outfile.write(str(pair)+"\n")

    return len(pixels[0])

def run_bright_pix(fitsfile, amps=image_utils.allAmps, outfile=None, \
                   verbose=True):
    """ Given an input FITS file, find bright pixels."""

    count_bright_pixels = 0
    for hdu in amps:
        pixels, columns, im = bright_pix(fitsfile, \
                                  image_utils.dm_hdu(hdu))
        count_bright_pixels += segment_results(hdu, pixels, outfile, verbose)

    ccdTotal = 'Total CCD Bright Pixels: ' + str(count_bright_pixels)
    if outfile != None:
        outfile.write(ccdTotal+"\n")

    if verbose:
        print ccdTotal


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
    count_bright_pixels = 0
    for hdu in image_utils.hdu_dict:
        pixels, columns, im = bright_pix(fitsfile, image_utils.dm_hdu(hdu))
        count_bright_pixels += segment_results(hdu, pixels)
    print "Total CCD Bright Pixels: ", count_bright_pixels

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Find bright pixels.')
    parser.add_argument('-i', '--infile', \
                        help="image file to be used for analysis", type=str)
    parser.add_argument('-a', '--amps', \
                        help="amps to be analyzed, separated by a space", \
                        type=int, nargs = '+', default=range(1, 17))
    parser.add_argument('-t', '--test', help="run test only", \
                        action='store_true', default=False)
    parser.add_argument('-o', '--outfile', \
                        help="output file to store results", type=str)
    parser.add_argument('-v', '--verbose', help="turn verbosity on", \
                        action='store_true', default=False)
    args = parser.parse_args()

    if args.outfile:
        outfile = open(args.outfile, 'w')
    else:
        outfile = None

    if args.test:
        run_test()
    else:
        run_bright_pix(args.infile, args.amps, outfile, args.verbose)

    if outfile:
        outfile.close()

