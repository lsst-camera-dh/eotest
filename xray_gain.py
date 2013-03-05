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
import image_utils as imUtils
from fe55_yield import fe55_yield

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
    def __init__(self, imfile, bbox=afwGeom.Box2I(), hdu=0, 
                 metadata=dafBase.PropertySet(), ccdtemp=None):
        self.image = afwImage.ImageF(imfile, hdu, metadata, bbox)
        self.fe55_yield = fe55_yield(ccdtemp)[0]
        #
        # Store numpy array of full segment (i.e., for an empty bbox,
        # which is the default in the ImageF constructor) for
        # self._footprint_signal, since footprint spans use full
        # segment image pixel coordinates.
        #
        self.arr = afwImage.ImageF(imfile, hdu).getArray()
        stats = afwMath.makeStatistics(self.image,
                                       afwMath.STDEVCLIP | afwMath.MEDIAN)
        self.noise = stats.getValue(afwMath.STDEVCLIP)
        self.median = stats.getValue(afwMath.MEDIAN)
        self.fp_sets = []
#    def _footprint_signal(self, footprint):
#        spans = footprint.getSpans()
#        total = 0
#        for span in spans:
#            total += sum(self.arr[span.getY()][span.getX0():span.getX1()+1])
#        return total - footprint.getNpix()*self.median
    def _footprint_signal(self, footprint, buff=1):
        bbox = footprint.getBBox()
        xmin = max(imUtils.imaging.getMinX(), bbox.getMinX() - buff)
        xmax = min(imUtils.imaging.getMaxX(), bbox.getMaxX() + buff)
        ymin = max(imUtils.imaging.getMinY(), bbox.getMinY() - buff)
        ymax = min(imUtils.imaging.getMaxY(), bbox.getMaxY() + buff)
        subarr = self.arr[ymin:ymax+1, xmin:xmax+1] 
        signal = sum(subarr.flat) - self.median*len(subarr.flat)
        return signal
    def gain(self, jmargin=None, max_npix=9, buff=1, regfile=None):
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
            signals = [self._footprint_signal(fp, buff) for fp in 
                       fp_set.getFootprints() if fp.getNpix() < max_npix]
            try:
                stats = afwMath.makeStatistics(signals, afwMath.MEANCLIP)
                mean = stats.getValue()
                values.append(self.fe55_yield/mean)
            except:
                pass
        my_gain = np.median(values)
        try:
            imed = np.where(values == my_gain)[0][0]
        except IndexError:
            values.pop()
            imed = np.where(values == np.median(values))[0][0]
        if regfile is not None:
            make_region_file(self.fp_sets[imed], regfile)
        return my_gain

def hdu_gains(infile, bbox=afwGeom.Box2I()):
    try:
        ccdtemp = afwImage.readMetadata(infile, 1).get('CCDTEMP')
    except:
        ccdtemp = -100.
    gains = {}
    for amp in imUtils.allAmps:
        fe55 = Fe55Gain(infile, bbox=bbox, hdu=imUtils.dm_hdu(amp),
                        ccdtemp=ccdtemp)
        gains[amp] = fe55.gain()
    return gains

if __name__ == '__main__':
    from simulation.sim_tools import SegmentExposure, writeFits

    ntrials = 10
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
