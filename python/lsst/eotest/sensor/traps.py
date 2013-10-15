"""
@brief Find traps from a pocket-pumped CCD frame.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import lsst.afw.detection as afwDetect
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.eotest.image_utils as imutils
from MaskedCCD import MaskedCCD
import pylab_plotter as plot

def traps(ccd, gains, outfile=None, cycles=100, nsig=6, amps=imutils.allAmps):
    if outfile is None:
        output = open(os.devnull, 'w')
    else:
        output = open(outfile, 'w')
    output.write('Amp  #x-pixel  y-pixel  trap_size\n')
    traps = {}
    for amp in amps:
        traps[amp] = []
        image = ccd.unbiased_and_trimmed_image(amp)
        stats = afwMath.makeStatistics(image,
                                       afwMath.MEDIAN | afwMath.STDEVCLIP,
                                       ccd.stat_ctrl)
        median = stats.getValue(afwMath.MEDIAN)
        stdev = stats.getValue(afwMath.STDEVCLIP)
        threshold = afwDetect.Threshold(median + nsig*stdev)
        fpset = afwDetect.FootprintSet(image, threshold)

        imaging = ccd.seg_regions[amp].imaging
        for fp in fpset.getFootprints():
            peak = [pk for pk in fp.getPeaks()][0]
            x, y = (peak.getIx() - imaging.getMinX(),
                    peak.getIy() - imaging.getMinY())
            bbox = afwGeom.Box2I(afwGeom.Point2I(x, imaging.getBeginY()),
                                 afwGeom.Extent2I(1, imaging.getHeight()))
            column = image.Factory(image, bbox)
            colmean = imutils.mean(column)
            trap_size = (peak.getPeakValue() - colmean)*gains[amp]/cycles
            traps[amp].append((peak.getIx(), peak.getIy(), trap_size))
            output.write(' %s  ' % imutils.channelIds[amp])
            output.write('%5i    %5i       %.1f\n' % traps[amp][-1])
            output.flush()
    output.close()
    return traps

if __name__ == '__main__':
    gains = dict([(amp, 5) for amp in imutils.allAmps])
    outfile = 'trap_list.txt'
    
    infile = 'work/sensorData/000-00/trap/debug/000-00_trap_ppump_debug.fits'
    mask_files = ()
    ccd = MaskedCCD(infile, mask_files=mask_files)

    my_traps = traps(ccd, gains, outfile=outfile)
