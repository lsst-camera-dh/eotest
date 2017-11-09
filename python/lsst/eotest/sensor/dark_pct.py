"""
@brief Function to compute the max dark current for a given percentile
of pixels.
"""
import numpy as np
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom

from image_utils import fits_median, unbias_and_trim


def dark_pct(files, percentile=90., hdu=2, gain=1):
    if percentile < 0 or percentile > 100:
        raise RuntimeError("percentile must be between 0 and 100")

    # Read exposure time from the primary HDU of the first file.
    exptime = afwImage.readMetadata(files[0], 1).get('EXPTIME')

    im = unbias_and_trim(fits_median(files))
    im *= gain
    im /= exptime

    npix = im.getHeight()*im.getWidth()
    imarr = np.sort(im.getArray().reshape(npix))

    return imarr[int(npix*float(percentile)/100.)]


if __name__ == '__main__':
    import glob
    files = glob.glob('data/dark*.fits')
    print dark_pct(files)
