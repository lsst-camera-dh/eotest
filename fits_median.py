"""
@brief Find the median image given a list of FITS filenames and hdu.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import numpy as num
import lsst.afw.image as afwImage

def fits_median(files, fix=1, hdu=2):
    ims = [afwImage.ImageF(f, hdu) for f in files]
    exptimes = [afwImage.readMetadata(f, 1).get('EXPTIME') for f in files]

    if min(exptimes) != max(exptimes):
        raise RuntimeError("Error: unequal exposure times")

    if fix == 1:
        medians = num.array([num.median(im.getArray()) for im in ims])
        med = sum(medians)/len(medians)
        errs = medians - med
        for im, err in zip(ims, errs):
            im -= err

    imcube = num.array([im.getArray() for im in ims])
    medim = afwImage.ImageF(num.median(imcube, axis=0))

    return medim

if __name__ == '__main__':
    import glob
    files = glob.glob('data/dark*.fits')
    im = fits_median(files)
