"""
@brief Find bright pixels and bright columns above a 
threshold = mean + nsig*sigma
"""
import numpy as np
import lsst.afw.image as afwImage

from image_utils import unbias_and_trim
from sim_tools import SegmentExposure, writeFits

def bright_pix(infile, hdu=2, nsig=5):
    im = unbias_and_trim(afwImage.ImageF(infile, hdu))
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

def write_test_image(outfile):
    seg = SegmentExposure()
    seg.add_bias(1e4, 10)
    seg.add_dark_current(300)
    seg.expose_flat(200)
    cols = seg.add_bright_cols(ncols=1, nsig=10)
    pix = seg.add_bright_pix(npix=100, nsig=10)
    writeFits((seg,), outfile)

if __name__ == '__main__':
    fitsfile = 'test_image.fits'
    write_test_image(fitsfile)
    pixels, columns, im = bright_pix(fitsfile)
