"""
@brief Module to perform standard operations on sensor images such
computing median images, unbiasing using the serial overscan region,
trimming, etc..
"""
import numpy as np
import numpy.random as random
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath

allAmps = range(1,17)

# Segment ID to HDU number in FITS dictionary
hdu_dict = dict( [ (1,'Segment10'), (2,'Segment11'), (3,'Segment12'),
                   (4,'Segment13'), (5,'Segment14'), (6,'Segment15'),
                   (7,'Segment16'), (8,'Segment17'), (9,'Segment07'),
                   (10,'Segment06'), (11,'Segment05'), (12,'Segment04'),
                   (13,'Segment03'), (14,'Segment02'), (15,'Segment01'),
                   (16,'Segment00') ] )

channelIds = dict([(i, hdu_dict[i][-2:]) for i in allAmps])

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

mean = lambda x : afwMath.makeStatistics(x, afwMath.MEAN).getValue()
median = lambda x : afwMath.makeStatistics(x, afwMath.MEDIAN).getValue()
stdev = lambda x : afwMath.makeStatistics(x, afwMath.STDEV).getValue()

detsize = '[0:4336,0:4044]'

def detsec(amp, dx=542, dy=2022):
    """DETSEC header keyword value for iraf mosaicking"""
    namps = len(allAmps)
    if amp < allAmps[namps/2]:
        x1 = dx*amp
        x2 = x1 - dx
        y1, y2 = 0, dy
    else:
        x1 = dx*(namps - amp)
        x2 = x1 + dx
        y1, y2 = 2*dy, dy
    return '[%i:%i,%i:%i]' % (x1, x2, y1, y2)

def dm_hdu(hdu):
    """ Compute DM HDU from the actual FITS file HDU."""
    return hdu+1

def bias(im, overscan=serial_overscan):
    """Compute the bias from the mean of the pixels in the serial
    overscan region."""
    return mean(im.Factory(im, overscan))

def bias_func(im, overscan=serial_overscan, fit_order=1):
    """Compute the bias by fitting a polynomial (linear, by default)
    to the mean of each row of the selected overscan region.  This
    returns a numpy.poly1d object that returns the fitted bias as
    function of pixel row."""
    try:
        imarr = im.Factory(im, overscan).getArray()
    except AttributeError: # Dealing with a MaskedImage
        imarr = im.Factory(im, overscan).getImage().getArray()
    ny, nx = imarr.shape
    rows = np.arange(ny)
    values = np.array([np.mean(imarr[j]) for j in rows])
    return np.poly1d(np.polyfit(rows, values, fit_order))

def bias_image(im, overscan=serial_overscan, fit_order=1):
    my_bias = bias_func(im, serial_overscan, fit_order)
    biasim = afwImage.ImageF(im.getDimensions())
    imarr = biasim.getArray()
    ny, nx = imarr.shape
    for row in range(ny):
        imarr[row] += my_bias(row)
    return biasim

def trim(im, imaging=imaging):
    "Trim the prescan and overscan regions."
    return im.Factory(im, imaging)

def unbias_and_trim(im, overscan=serial_overscan, imaging=imaging,
                    apply_trim=True, fit_order=1):
    """Subtract bias calculated from overscan region and optionally trim 
    prescan and overscan regions."""
#    imarr = im.getArray()
#    ny, nx = imarr.shape
#
#    my_bias = bias_func(im, overscan, fit_order)
#    for row in range(ny):
#        imarr[row] -= my_bias(row)
    im -= bias_image(im, serial_overscan, fit_order)
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

class SubRegionSampler(object):
    def __init__(self, dx, dy, nsamp, imaging):
        self.dx = dx
        self.dy = dy
        #
        # Generate sub-regions at random locations on the segement
        # imaging region.
        #
        self.xarr = random.randint(imaging.getWidth() - dx - 1, size=nsamp)
        self.yarr = random.randint(imaging.getHeight() - dy - 1, size=nsamp)
    def bbox(self, x, y):
        return afwGeom.Box2I(afwGeom.Point2I(int(x), int(y)),
                             afwGeom.Extent2I(self.dx, self.dy))
    def subim(self, im, x, y):
        return im.Factory(im, self.bbox(x, y))

if __name__ == '__main__':
    import glob

    files = glob.glob('data/dark*.fits')
    im = fits_median(files)

    bar = unbias_and_trim(im)
