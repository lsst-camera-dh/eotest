import numpy as num
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath

infile = 'fe55_0060s_000.fits'
pct = 90.
hdu = 2
threshold = 1000

im = afwImage.ImageF(infile, hdu)

bbox = afwGeom.Box2I(afwGeom.Point2I(10, 0), afwGeom.Point2I(521, 2002))
overscan = afwGeom.Box2I(afwGeom.Point2I(522, 0), 
                         afwGeom.Point2I(im.getWidth()-1, im.getHeight()-1))
bias = num.mean(im.Factory(im, overscan).getArray())

imarr = im.getArray()

col_means = [num.mean(imarr[i,:]) for i in range(512)]

bright_pix = num.where(
