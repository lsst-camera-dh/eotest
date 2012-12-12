"""
@brief Decorator class to abstract sensor properties such as 
overscan, bias, active region, etc.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import numpy as num
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage

class SensorImage(object):
    overscan = afwGeom.Box2I(afwGeom.Point2I(525, 100), 
                             afwGeom.Extent2I(5, 1900))
    active = afwGeom.Box2I(afwGeom.Point2I(10, 0), 
                           afwGeom.Point2I(521, 2000))
    def __init__(self, im, trim=False):
        self._im = im
        self.bias = num.mean(im.Factory(im, self.overscan).getArray())
        self.npix = im.getHeight()*im.getWidth()
        if trim:
            self.image = im.Factory(im, self.active)
        else:
            self.image = im
    def unbias(self):
        self.image -= self.bias
        return self
    def Factory(self, sensor_image, bbox):
        return self.image.Factory(sensor_image.image, bbox)
    def __getattr__(self, attr):
        # Delegate to Image object.
        return getattr(self.image, attr)

if __name__ == '__main__':
    im = afwImage.ImageF('fe55_0060s_000.fits')

    foo = SensorImage(im)
