"""
@brief Module to perform standard operations on sensor images such
computing median images, unbiasing using the overscan region,
trimming, etc..
"""
import numpy as np
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage

overscan = afwGeom.Box2I(afwGeom.Point2I(525, 100), 
                         afwGeom.Extent2I(5, 1900))
active = afwGeom.Box2I(afwGeom.Point2I(10, 0), 
                       afwGeom.Point2I(521, 2000))

def fits_median(files, hdu=2, fix=True):
    """Compute the median image from a set of image FITS files."""
    ims = [afwImage.ImageF(f, hdu) for f in files]
    exptimes = [afwImage.readMetadata(f, 1).get('EXPTIME') for f in files]

    if min(exptimes) != max(exptimes):
        raise RuntimeError("Error: unequal exposure times")

    if fix:
        medians = np.array([np.median(im.getArray()) for im in ims])
        med = sum(medians)/len(medians)
        errs = medians - med
        for im, err in zip(ims, errs):
            im -= err

    imcube = np.array([im.getArray() for im in ims])
    medim = afwImage.ImageF(np.median(imcube, axis=0))

    return medim

def unbias_and_trim(im, active=active, overscan=overscan, trim=True):
    """Subtract bias calculated from overscan region and optionally trim 
    prescan and overscan regions."""
    bias = np.mean(im.Factory(im, overscan).getArray())
    im -= bias
    if trim:
        return im.Factory(im, active)
    else:
        return im

if __name__ == '__main__':
    import glob

    files = glob.glob('data/dark*.fits')
    im = fits_median(files)

    bar = unbias_and_trim(im)
