import numpy as np
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath

from image_utils import unbias_and_trim

infile = 'fe55_0060s_000.fits'
pct = 90.
hdu = 2
threshold = 1000

im = unbias_and_trim(afwImage.ImageF(infile, hdu))

imarr = im.getArray()

# Mean value of each CCD column.
col_means = [np.mean(imarr[i,:]) for i in range(512)]

# Find columns where mean > pct/100*(image mean)
bcols = np.where(col_means >= (mean(imarr)*pct/100.))
count = len(bcols[0])

bright_pix = np.where(imarr > threshold)
bp_count = len(bright_pix[0])

if bp_count > 0:
    bp_mean = np.mean(imarr[bright_pix])
else:
    bp_mean = 0

# Find bright columns
colavgs = [np.median(col) for col in imarr]
