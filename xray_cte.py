import glob
import lsst.afw.detection as afwDetect
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import image_utils as imUtils
import numpy as np
import numpy.linalg as linalg
try:
    import pylab_plotter as plot
except ImportError:
    plot = None

class XrayCte(object):
    def __init__(self, image):
        self.image = image
        self.image -= imUtils.bias(image)
        self.imarr = image.getArray()

        stat_control = afwMath.MEANCLIP | afwMath.STDEVCLIP 
        stats = afwMath.makeStatistics(image, stat_control)
        self.mean = stats.getValue(afwMath.MEANCLIP)
        self.stdev = stats.getValue(afwMath.STDEVCLIP)
    def find_hits(self, nsig=2, DN_range=None, make_plots=True):
        threshold = afwDetect.Threshold(self.mean + 2*self.stdev)
        fpset = afwDetect.FootprintSet(self.image, threshold)
        footprints = [fp for fp in fpset.getFootprints()]
        zarr, xarr, yarr = [], [], []
        for fp in footprints:
            if fp.getNpix() < 9:
                f, ix, iy = self._footprint_info(fp)
                if DN_range is None or DN_range[0] < f < DN_range[1]:
                    zarr.append(f)
                    xarr.append(ix)
                    yarr.append(iy)
        self.xarr = np.array(xarr)
        self.yarr = np.array(yarr)
        self.zarr = np.array(zarr)
        if plot is not None and make_plots:
            if DN_range is not None:
                sigrange = DN_range
            else:
                sigrange = (0, 1000)
            plot.histogram(zarr, xrange=sigrange,
                           xname='DN', yname='counts / bin')
            plot.xyplot(xarr, zarr, yrange=sigrange,
                        xname='x pixel', yname='DN')
    def _footprint_info(self, fp):
        imarr = self.imarr
        spans = fp.getSpans()
        total = 0
        for span in spans:
            total += sum(imarr[span.getY()][span.getX0():span.getX1()+1])
        peakvalues = np.array([x.getPeakValue() for x in fp.getPeaks()])
        peaks = [x for x in fp.getPeaks()]
        if len(peakvalues) > 1:
            indx = np.where(peakvalues == max(peakvalues))
            ii = indx[0][0]
        else:
            ii = 0
        ix, iy = peaks[ii].getIx(), peaks[ii].getIy()
        return total - self.mean*fp.getNpix(), ix, iy
    def fit_1d(self, xmin=100, make_plots=True):
        # Omit pixels affected by edge roll-off.
        indx = np.where(self.xarr > xmin)
        xarr = self.xarr[indx]
        yarr = self.yarr[indx]
        zarr = self.zarr[indx]

        ax, bx = np.polyfit(xarr, zarr, 1)
        CTIx = -ax/bx
        print "1D CTI estimates:"
        print "x-direction:"
        print "   CTI =", CTIx
        print "   gain =", 1620./bx
        if plot is not None and make_plots:
            fx = np.poly1d((ax, bx))
            xx = np.linspace(0, 600, 100)
            plot.xyplot(xarr, zarr, yrange=(300, 400),
                        xname='x pixel', yname='DN')
            plot.curve(xx, fx(xx), oplot=1, lineStyle='--', color='r')

        ay, by = np.polyfit(yarr, zarr, 1)
        CTIy = -ay/by
        print "y-direction:"
        print "   CTI =", CTIy
        print "   gain =", 1620./by
        if plot is not None and make_plots:
            fy = np.poly1d((ay, by))
            yy = np.linspace(0, 2500, 100)
            plot.xyplot(yarr, zarr, yrange=(300, 400),
                        xname='y pixel', yname='DN')
            plot.curve(yy, fy(yy), oplot=1, lineStyle='--', color='r')
        print
    def fit_2d(self, xmin=100):
        # Omit pixels affected by edge roll-off.
        indx = np.where(self.xarr > xmin)
        xarr = self.xarr[indx]
        yarr = self.yarr[indx]
        zarr = self.zarr[indx]
        
        A = np.matrix([[sum(xarr**2), sum(xarr*yarr), -sum(xarr)],
                       [sum(xarr*yarr), sum(yarr**2), -sum(yarr)],
                       [-sum(xarr), -sum(yarr), len(xarr)]])
        B = [-sum(xarr*zarr), -sum(yarr*zarr), sum(zarr)]
    
        CTIx, CTIy, sig0 = linalg.solve(A, B)
        print "2D CTI estimates:"
        print "   CTI(x) =", CTIx/sig0
        print "   CTI(y) =", CTIy/sig0
        print "   gain =", 1620./sig0

if __name__ == '__main__':
    fe55 = '/nfs/farm/g/lsst/u1/testData/eotestData/000-00/xray/data/000_00_fe55_0600s_008.fits'
    for amp in imUtils.allAmps:
        image = afwImage.ImageF(fe55, imUtils.dm_hdu(amp))
        print "Segment %s" % imUtils.channelIds[amp]
        cte = XrayCte(image)
        cte.find_hits(nsig=2, DN_range=(330, 360))
        cte.fit_1d(200)
        cte.fit_2d(200)
