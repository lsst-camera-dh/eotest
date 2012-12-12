import glob
import numpy as num
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom

from fits_median import fits_median
from sensorImage import SensorImage

def dark_pct(files, percentile=90., hdu=2, gain=1):
    if percentile < 0 or percentile > 100:
        raise RuntimeError("percentile must be between 0 and 100")

    med_image = fits_median(files)
    exptime = afwImage.readMetadata(files[0], 1).get('EXPTIME')
    im = SensorImage(med_image).unbias()
    im.image *= gain
    im.image /= exptime
    npix = im.getHeight()*im.getWidth()
    imarr = num.sort(im.getArray().reshape(npix))
    return imarr[int(npix*float(percentile)/100.)]

if __name__ == '__main__':
    files = glob.glob('data/dark*.fits')
    print dark_pct(files)

