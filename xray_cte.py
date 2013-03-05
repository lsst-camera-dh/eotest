"""
@brief Compute charge transfer inefficiency (CTI) by fitting Fe55 data
using linear functions, either separately in 1D, for parallel and
serial directions, or for both directions at once in a 2D fit.  These
fits also give estimates of the system gain.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import sys
import glob
import numpy as np
import numpy.linalg as linalg
import lsst.afw.detection as afwDetect
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import image_utils as imUtils
from fe55_yield import fe55_yield
try:
    import pylab_plotter as plot
except ImportError:
    plot = None

class XrayCte(object):
    def __init__(self, image, ccdtemp=None):
        if (image.getWidth() == imUtils.full_segment.getWidth() and
            image.getHeight() == imUtils.full_segment.getHeight()):
            # Trim to imaging area.
            self.image = image.Factory(image, imUtils.imaging)
        else:
            message = "XrayCte constructor: must pass an untrimmed image"
            raise RuntimeError(message)
        self.imarr = image.getArray()
        stat_control = afwMath.MEANCLIP | afwMath.STDEVCLIP 
        stats = afwMath.makeStatistics(image, stat_control)
        self.mean = stats.getValue(afwMath.MEANCLIP)
        self.stdev = stats.getValue(afwMath.STDEVCLIP)
        self.fe55_yield = fe55_yield(ccdtemp)[0]
    def find_hits(self, nsig=2, gain_range=(2, 10), buff=1, make_plots=True):
        DN_range = (self.fe55_yield/gain_range[1],
                    self.fe55_yield/gain_range[0])
        threshold = afwDetect.Threshold(self.mean + nsig*self.stdev)
        fpset = afwDetect.FootprintSet(self.image, threshold)
        zarr, xarr, yarr = [], [], []
        for fp in fpset.getFootprints():
            if fp.getNpix() < 9:
                f, ix, iy = self._footprint_info(fp, buff)
                if DN_range[0] < f < DN_range[1]:
                    zarr.append(f)
                    xarr.append(ix)
                    yarr.append(iy)
        self.xarr = np.array(xarr)
        self.yarr = np.array(yarr)
        self.zarr = np.array(zarr)
        median_signal = imUtils.median(self.zarr)
        self.median_signal = median_signal
        self.sigrange = (median_signal*(1. - 2./np.sqrt(self.fe55_yield)),
                         median_signal*(1. + 2./np.sqrt(self.fe55_yield)))
        self.sig5range = (median_signal*(1. - 5./np.sqrt(self.fe55_yield)),
                          median_signal*(1. + 5./np.sqrt(self.fe55_yield)))
        if plot is not None and make_plots:
            plot.histogram(self.zarr, xrange=self.sig5range,
                           yrange=(0, 200), bins=100,
                           xname='DN', yname='entries / bin')
            plot.vline(median_signal)
            plot.vline(self.sigrange[0], color='g')
            plot.vline(self.sigrange[1], color='g')
            plot.xyplot(xarr, zarr, yrange=DN_range,
                        xname='x pixel', yname='DN')
            plot.hline(self.sigrange[0], color='g')
            plot.hline(self.sigrange[1], color='g')
    def _footprint_info(self, fp, buff=1):
        bbox = fp.getBBox()
        xmin = max(imUtils.imaging.getMinX(), bbox.getMinX() - buff)
        xmax = min(imUtils.imaging.getMaxX(), bbox.getMaxX() + buff)
        ymin = max(imUtils.imaging.getMinY(), bbox.getMinY() - buff)
        ymax = min(imUtils.imaging.getMaxY(), bbox.getMaxY() + buff)
        subarr = self.imarr[ymin:ymax+1, xmin:xmax+1] 
        signal = sum(subarr.flat) - self.mean*len(subarr.flat)
        peakvalues = np.array([x.getPeakValue() for x in fp.getPeaks()])
        peaks = [x for x in fp.getPeaks()]
        if len(peakvalues) > 1:
            indx = np.where(peakvalues == max(peakvalues))
            ii = indx[0][0]
        else:
            ii = 0
        ix, iy = peaks[ii].getIx(), peaks[ii].getIy()
        return signal, ix, iy
    def fit_1d(self, xmin=100, make_plots=True):
        # Omit pixels affected by edge roll-off and that are outside
        # the nominal signal range.
        indx = np.where((self.xarr > xmin) & (self.zarr > self.sigrange[0]) &
                        (self.zarr < self.sigrange[1]))
        xarr = self.xarr[indx]
        yarr = self.yarr[indx]
        zarr = self.zarr[indx]

        ax, bx = np.polyfit(xarr, zarr, 1)
        CTIx = -ax/bx

        ay, by = np.polyfit(yarr, zarr, 1)
        CTIy = -ay/by

        if plot is not None and make_plots:
            self._plot1d(ax, bx, 'x pixel', xarr, zarr)
            self._plot1d(ay, by, 'y pixel', yarr, zarr)

        return CTIx, CTIy, self.fe55_yield/bx, self.fe55_yield/by
    def _plot1d(self, a, b, xlabel, xarr, zarr):
        f = np.poly1d((a, b))
        xx = np.linspace(0, 2500, 5)
        win = plot.xyplot(xarr, zarr, yrange=self.sig5range,
                          xname=xlabel, yname='DN')
        plot.curve(xx, f(xx), oplot=1, lineStyle='--', color='r')
        plot.hline(self.sigrange[0], color='g')
        plot.hline(self.sigrange[1], color='g')
        return win
    def fit_2d(self, xmin=100):
        # Omit pixels affected by edge roll-off and that are outside
        # the nominal signal range.
        indx = np.where((self.xarr > xmin) & (self.zarr > self.sigrange[0]) &
                        (self.zarr < self.sigrange[1]))
        xarr = self.xarr[indx]
        yarr = self.yarr[indx]
        zarr = self.zarr[indx]
        
        A = np.matrix([[sum(xarr**2), sum(xarr*yarr), -sum(xarr)],
                       [sum(xarr*yarr), sum(yarr**2), -sum(yarr)],
                       [-sum(xarr), -sum(yarr), len(xarr)]])
        B = [-sum(xarr*zarr), -sum(yarr*zarr), sum(zarr)]
    
        CTIx, CTIy, sig0 = linalg.solve(A, B)
        return CTIx/sig0, CTIy/sig0, self.fe55_yield/sig0

if __name__ == '__main__':
    fe55 = '/nfs/farm/g/lsst/u1/testData/eotestData/000-00/xray/data/000_00_fe55_0600s_000.fits'
#    fe55 = '/nfs/farm/g/lsst/u1/testData/SIMData/000-00/Fe55/Fe55_exp_000-00_00.fits'
    try:
        ccdtemp = afwImage.readMetadata(fe55, 1).get('CCDTEMP')
    except:
        ccdtemp = -100
    print "CCD temp:", ccdtemp
    print "Segment   serial CTI   parallel CTI   gain ests."
    make_plots = False
    for amp in imUtils.allAmps:
#    make_plots = True
#    for amp in (1,):
        image = afwImage.ImageF(fe55, imUtils.dm_hdu(amp))
        cte = XrayCte(image, ccdtemp=ccdtemp)
        cte.find_hits(nsig=5, make_plots=make_plots)
        sys.stdout.write("%s       " % imUtils.channelIds[amp])
        results = cte.fit_1d(200, make_plots=make_plots)
        print "%11.4e   %11.4e    %.3f  %.3f"  % results
        results = cte.fit_2d(200)
        print "         %11.4e   %11.4e    %.3f" % results
        
