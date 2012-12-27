"""
@brief Re-implementation of P. Doherty's IDL function to compute
sensor gain from Fe55 data.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import numpy as np
import lsst.afw.math as afwMath
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.detection as afwDetection
import lsst.daf.base as dafBase

_ds9_header = """# Region file format: DS9 version 4.0
global color=green font="helvetica 10 normal" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source
linear
"""

def make_region_file(fpset, outfile='ds9.reg'):
    output = open(outfile, 'w')
    output.write(_ds9_header)
    for fp in fpset.getFootprints():
        peak = fp.getPeaks()[0]
        output.write('point(%i,%i) # point=circle\n' 
                     % (peak.getIx(), peak.getIy()))
    output.close()

class Fe55Gain(object):
    _stat_controls = (afwMath.NPOINT | afwMath.MEAN | afwMath.STDEV 
                      | afwMath.MEDIAN)
    def __init__(self, imfile, bbox=afwGeom.Box2I(), hdu=0, 
                 metadata=dafBase.PropertySet()):
        self.image = afwImage.ImageF(imfile, hdu, metadata, bbox)
        self.arr = self.image.getArray()
        my_arr = self.arr.reshape(self.image.getHeight()*self.image.getWidth())
        stats = self._clipped_stats(my_arr.tolist())
        self.noise = stats.getValue(afwMath.STDEV)
        self.median = stats.getValue(afwMath.MEDIAN)
        self.fp_sets = []
    def _footprint_signal(self, footprint):
        spans = footprint.getSpans()
        total = 0
        for span in spans:
            total += sum(self.arr[span.getY()][span.getX0():span.getX1()+1])
        return total - footprint.getNpix()*self.median
    def _clipped_stats(self, xvals, nsig=3):
        stats = afwMath.makeStatistics(xvals, self._stat_controls)
        med = stats.getValue(afwMath.MEDIAN)
        stdev = stats.getValue(afwMath.STDEV)
        xmin, xmax = med - stdev*nsig, med + stdev*nsig
        my_xvals = [x for x in xvals if xmin <= x <= xmax]
        return afwMath.makeStatistics(my_xvals, self._stat_controls)
    def gain(self, jmargin=None, max_npix=9, regfile=None):
        if jmargin is None:
            jmargins = range(10)
        else:
            jmargins = (jmargin,)
        values = []
        for j in jmargins:
            margin = self.noise*(j + 3)
            threshold = afwDetection.Threshold(self.median + margin)
            fp_set = afwDetection.FootprintSet(self.image, threshold)
            self.fp_sets.append(fp_set)
            signals = [self._footprint_signal(fp) for fp in 
                       fp_set.getFootprints() if fp.getNpix() < max_npix]
            mean = self._clipped_stats(signals).getValue(afwMath.MEAN)
            values.append(1620./mean)
        imed = np.argsort(values)[len(values)/2]
        if regfile is not None:
            make_region_file(self.fp_sets[imed], regfile)
        return values[imed]

if __name__ == '__main__':
    from sim_tools import SegmentExposure, writeFits

    ntrials = 100
    gains = []
    hdu = 1
    for i in range(ntrials):
        seg = SegmentExposure(exptime=1, gain=3)
        seg.add_bias(level=2000, sigma=2)
        seg.add_dark_current(level=1000)
        seg.add_Fe55_hits(nxrays=200)
        imfile = writeFits((seg,), 'test_fe55_image.fits')

        f55 = Fe55Gain(imfile, hdu=hdu+1)
        gains.append(f55.gain())
        print i, gains[-1]
    print np.mean(gains), np.median(gains), np.std(gains)
