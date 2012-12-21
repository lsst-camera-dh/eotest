"""
@brief Module to perform standard operations on sensor images such
computing median images, unbiasing using the serial overscan region,
trimming, etc..
"""
import numpy as np
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage

full_segment = afwGeom.Box2I(afwGeom.Point2I(0, 0),
                             afwGeom.Point2I(541, 2021))

prescan = afwGeom.Box2I(afwGeom.Point2I(0, 0),
                        afwGeom.Point2I(9, 2021))

imaging = afwGeom.Box2I(afwGeom.Point2I(10, 0),
                        afwGeom.Point2I(521, 2001))

serial_overscan = afwGeom.Box2I(afwGeom.Point2I(522, 0), 
                                afwGeom.Point2I(541, 2021))

parallel_overscan = afwGeom.Box2I(afwGeom.Point2I(10, 2002),
                                  afwGeom.Point2I(522, 2021))

overscan = serial_overscan  # for backwards compatibility

def bias(im, overscan=serial_overscan):
    "Compute the bias from the serial overscan region."
    return np.mean(im.Factory(im, overscan).getArray())

def trim(im, imaging=imaging):
    "Trim the prescan and overscan regions."
    return im.Factory(im, imaging)

def unbias_and_trim(im, overscan=serial_overscan, imaging=imaging,
                    apply_trim=True):
    """Subtract bias calculated from overscan region and optionally trim 
    prescan and overscan regions."""
    im -= bias(im, overscan)
    if apply_trim:
        return trim(im, imaging)
    else:
        return im

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

if __name__ == '__main__':
    import glob

    files = glob.glob('data/dark*.fits')
    im = fits_median(files)

    bar = unbias_and_trim(im)
