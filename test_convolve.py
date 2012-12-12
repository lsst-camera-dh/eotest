"""
@brief Example use of convolution tool on a FITS image "directly",
i.e., avoiding the overhead associated with retrieving data with the
data butler.
"""
#
# $Header: /usr/local/CVS/SLAC/users/jchiang/lsst_dm_tools/test_convolve.py,v 1.3 2012/11/08 19:26:47 jchiang Exp $
#
import numpy as num
import lsst.afw.math as afwMath
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.detection as afwDetection
import lsst.afw.display.ds9 as ds9

imfile = 'imsim_886894611_R23_S11_C00_E000.fits.gz'

input_image = afwImage.MaskedImageF(imfile)

#
# One can create a delta function kernel directly.
#
deltaKernel = afwMath.DeltaFunctionKernel(1, 1, afwGeom.Point2I(0, 0))

#
# Using the PSF factory function seems like the best way to create one
# with an analytic function.
#
fwhm = 5
psf = afwDetection.createPsf("DoubleGaussian", 15, 15, 
                             fwhm/(2*num.sqrt(2*num.log(2))))
gaussKernel = psf.getKernel()

#
# afwMath.convolve is swig-exposed C++, so it has the unpythonic
# interface wherein the output image is passed as an argument.
#
output_image = input_image.Factory(input_image.getDimensions())
afwMath.convolve(output_image, input_image, gaussKernel)

#
# Display unconvolved and convolved images with ds9 and blink.
#
ds9.mtv(input_image)
ds9.incrDefaultFrame()
ds9.mtv(output_image)
ds9.ds9Cmd('blink')

