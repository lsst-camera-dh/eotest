"""
@brief Find traps from a pocket-pumped CCD frame.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import numpy as np
import pyfits
import lsst.afw.detection as afwDetect
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.eotest.image_utils as imutils
from MaskedCCD import MaskedCCD
import pylab_plotter as plot

def traps(ccd, gains, outfile='traps.fits', cycles=100, nsig=6, 
          amps=imutils.allAmps):
    """
    Find traps from a pocket-pumped image and write the results to a
    FITS file as a binary table.
    """
    traps = {}
    row = 0
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
    # Write the results as a FITS binary table.
    nrows = sum([len(traps[amp]) for amp in traps])
    output = pyfits.HDUList()
    output.append(pyfits.PrimaryHDU())
    colnames = ['AMPLIFIER', 'XPOS', 'YPOS', 'TRAP_SIZE']
    formats = 'IIII'
    units = ['None', 'pixel', 'pixel', 'electrons']
    columns = [np.zeros(nrows, dtype=int)]*4
    hdu = pyfits.new_table([pyfits.Column(name=colname, format=format,
                                          unit=unit, array=column) 
                            for colname, format, unit, column in
                            zip(colnames, formats, units, columns)])
    hdu.name = 'TRAPS'
    output.append(hdu)
    row = 0
    for amp in amps:
        for xpos, ypos, trap_size in traps[amp]:
            output['TRAPS'].data[row]['AMPLIFIER'] = amp
            output['TRAPS'].data[row]['XPOS'] = xpos
            output['TRAPS'].data[row]['YPOS'] = ypos
            output['TRAPS'].data[row]['TRAP_SIZE'] = trap_size
            row += 1
    output.writeto(outfile, clobber=True)

    return traps

if __name__ == '__main__':
    gains = dict([(amp, 5) for amp in imutils.allAmps])
    outfile = 'trap_list.fits'
    
    infile = 'work/sensorData/000-00/trap/debug/000-00_trap_ppump_debug.fits'
    mask_files = ()
    ccd = MaskedCCD(infile, mask_files=mask_files)

    my_traps = traps(ccd, gains, outfile=outfile)
