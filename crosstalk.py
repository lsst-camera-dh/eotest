"""
@brief Crosstalk calculations for system and detector inputs.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import time
import numpy as np
import pylab
import lsst.afw.detection as afwDetect
import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath
import image_utils as imutils
from MaskedCCD import MaskedCCD
from BrightPixels import BrightPixels

def get_stats(raw_image, stat_ctrl, imaging=imutils.imaging):
    image = imutils.unbias_and_trim(raw_image, imaging=imaging)
    flags = afwMath.MEDIAN | afwMath.STDEV
    stats = afwMath.makeStatistics(image, flags, stat_ctrl)
    return stats.getValue(afwMath.MEDIAN), stats.getValue(afwMath.STDEV)

def aggressor(ccd):
    """Guess the aggressor amp based on the maximum pixel value."""
    max_pix = lambda amp : max(ccd[amp].getImage().getArray().flat)
    candidate = imutils.allAmps[0]
    max_pix_val = max_pix(candidate)
    for amp in imutils.allAmps[1:]:
        val = max_pix(amp)
        if val > max_pix_val:
            candidate = amp
            max_pix_val = val
    return candidate, max_pix_val

def column_mean(raw_image, col, stat_ctrl, imaging=imutils.imaging):
    reg = afwGeom.Box2I(afwGeom.Point2I(col, imaging.getMinY()),
                        afwGeom.Extent2I(1, imaging.getHeight()))
    image = imutils.unbias_and_trim(raw_image, imaging=imaging)
    subim = image.Factory(image, reg)
    flags = afwMath.MEAN | afwMath.STDEV | afwMath.NPOINT
    stats = afwMath.makeStatistics(subim, flags, stat_ctrl) 
    mean = stats.getValue(afwMath.MEAN)
    npts = stats.getValue(afwMath.NPOINT)
    sigma = stats.getValue(afwMath.STDEV)/np.sqrt(npts)
    return np.array((mean, sigma))

def system_crosstalk(ccd, aggressor_amp, dnthresh=None, nsig=5):
    """
    Compute the system crosstalk.  dnthresh is the threshold in DN for
    detecting the illuminated column in the aggressor amplifier; if
    set to None, then nsig*clipped_stdev above median is used for
    the threshold.
    
    This routine is optimized for system crosstalk input, i.e., for
    which the aggressor amp has a single illuminated column.
    """
    #
    # Find bright column of aggressor amplifier.
    #
    # Set exptime and gain to unity so that the e- threshold used by
    # BrightPixels converts directly to DN.
    #
    exptime = 1
    gain = 1
    if dnthresh is None:
        median, stdev = get_stats(ccd[aggressor_amp], ccd.stat_ctrl)
        dnthresh = median + nsig*stdev
    bp = BrightPixels(ccd[aggressor_amp], exptime=exptime, gain=gain,
                      ethresh=dnthresh)
    pixels, columns = bp.find()

    if len(columns) > 1:
        raise RuntimeError("More than one aggressor column found.")

    agg_col = columns[0]
    agg_mean = column_mean(ccd[aggressor_amp], agg_col, ccd.stat_ctrl)[0]
    
    ratios = dict([(amp, column_mean(ccd[amp], agg_col, ccd.stat_ctrl)/agg_mean)
                   for amp in imutils.allAmps])
    return ratios

def get_footprint(fp_set):
    footprints = [fp for fp in fp_set.getFootprints()]
    if len(footprints) > 1:
        message = "More than one spot image found in aggressor amplifier:\n"
        for i, fp in enumerate(footprints):
            message += '%3i  %6i\n' % (i, [x.getPeakValue() 
                                           for x in fp.getPeaks()][0])
        raise RuntimeError(message)
    fp = footprints[0]
    peak_value = max([x.getPeakValue() for x in fp.getPeaks()])
    return fp, peak_value

def extract_mean_signal_2(masked_image, footprint, stat_ctrl):
    masked_image -= imutils.bias_image(masked_image)
    stdev = afwMath.makeStatistics(masked_image, afwMath.STDEVCLIP).getValue()
    signal = 0
    npix = 0
    for span in footprint.getSpans():
        width = span.getX1() - span.getX0() + 1
        bbox = afwGeom.Box2I(afwGeom.Point2I(span.getX0(), span.getY()),
                             afwGeom.Extent2I(width, 1))
        subim = masked_image.Factory(masked_image, bbox)
        stats = afwMath.makeStatistics(subim, afwMath.SUM | afwMath.NPOINT,
                                       stat_ctrl)
        signal += stats.getValue(afwMath.SUM)
        npix += stats.getValue(afwMath.NPOINT)
    return np.array((signal/npix, stdev))

def extract_mean_signal(masked_image, footprint, stat_ctrl):
    image = masked_image.getImage()
    maskarr = masked_image.getMask().getArray()
    image -= imutils.bias_image(image)
    stdev = afwMath.makeStatistics(image, afwMath.STDEVCLIP).getValue()
    imarr = image.getArray()
    signal = 0
    npix = 0
    for span in footprint.getSpans():
        y = span.getY()
        for x in range(span.getX0(), span.getX1()+1):
            npix += 1
            if maskarr[y][x] == 0:
                signal += imarr[y][x]
    return np.array((signal/float(npix), stdev))

def detector_crosstalk(ccd, aggressor_amp, dnthresh=None, nsig=5,
                       signal_extractor=extract_mean_signal):
    """
    Compute detector crosstalk from a spot image in the aggressor
    amplifier. dnthresh is the threshold in DN for detecting the
    illuminated column in the aggressor amplifier; if set to None,
    then nsig*clipped_stdev above median is used for the threshold.
    """
    image = imutils.unbias_and_trim(ccd[aggressor_amp])
    #
    # Extract footprint of spot image using nominal detection
    # threshold.
    #
    if dnthresh is None:
        median, stdev = get_stats(ccd[aggressor_amp], ccd.stat_ctrl)
        dnthresh = median + nsig*stdev
    threshold = afwDetect.Threshold(dnthresh)
    fp_set = afwDetect.FootprintSet(image, threshold)
    footprint, peak_value = get_footprint(fp_set)

    agg_mean = signal_extractor(ccd[aggressor_amp], footprint, ccd.stat_ctrl)[0]
    ratios = dict([(amp, signal_extractor(ccd[amp], footprint, ccd.stat_ctrl)
                    /agg_mean) for amp in imutils.allAmps])
    return ratios

class CrosstalkMatrix(object):
    def __init__(self, filename=None, namps=16):
        self.filename = filename
        self.namps = namps
        if self.filename is not None and os.path.isfile(self.filename):
            self._read_matrix()
        else:
            self._set_matrix()
    def set_row(self, agg, ratios):
        self.matrix[agg-1] = np.array([ratios[amp][0] for amp 
                                       in imutils.allAmps])
    def _set_matrix(self):
        self.matrix = np.zeros((self.namps, self.namps), dtype=np.float)
    def _read_matrix(self):
        self._set_matrix()
        input = open(self.filename, 'r')
        amp = 0
        for line in input:
            if line[0] == '#':
                continue
            self.matrix[amp] = np.array([float(x) for x 
                                         in line.strip().split()])
            amp += 1
    def write(self, outfile=None):
        if outfile is None:
            outfile = self.filename
        else:
            self.filename = outfile
        output = open(outfile, 'w')
        output.write('#')
        for amp in range(1, self.namps+1, 1):
            output.write('%02i  ' % amp)
        output.write('\n')
        for agg in range(self.namps):
            for victim in range(self.namps):
                output.write('%12.4e  ' % self.matrix[agg][victim])
            output.write('\n')
        output.close()
    def plot_matrix(self, title=None, cmap_range=(0.6, 0.4), precision=3, 
                    scale_factor=1e2, fontsize=10, figsize=(10, 6)):
        cmin, cmax = cmap_range
        my_matrix = np.copy(self.matrix)*scale_factor
        for i in range(self.namps):
            my_matrix[i][i] = 0
        thresh = 10**(-precision)
        cdict = dict(red=((0, cmin, cmin), (1, cmax, cmax)),
                     green=((0, cmin, cmin), (1, cmax, cmax)),
                     blue=((0, cmin, cmin), (1, cmax, cmax)))
        unsat_grey = pylab.matplotlib.colors.LinearSegmentedColormap('my_grey', 
                                                                     cdict, 256)
        fig = pylab.figure(figsize=figsize)
        axes = fig.add_subplot(111)
        pylab.imshow(np.abs(my_matrix), interpolation='nearest',
                     cmap=unsat_grey, aspect='auto')
        axes.set_xticks(range(self.namps))
        axes.set_yticks(range(self.namps))
        axes.set_xticklabels(['%i' % i for i in range(1, self.namps+1)])
        axes.set_yticklabels(['%i' % i for i in range(1, self.namps+1)])
        if title is not None:
            axes.set_title(title)
        ny, nx = my_matrix.shape
        for x in range(nx):
            for y in range(ny):
                if np.abs(my_matrix[y][x]) < thresh:
                    label = '0'
                else:
                    label = '%.3f' % my_matrix[y][x]
                pylab.text(x, y, label, horizontalalignment='center', 
                           fontsize=fontsize)
    def __sub__(self, other):
        result = CrosstalkMatrix()
        result.matrix = self.matrix - other.matrix
        return result
    def __add__(self, other):
        result = CrosstalkMatrix()
        result.matrix = self.matrix + other.matrix
        return result

def make_crosstalk_matrix(file_list, mask_files=(),
                          extractor=detector_crosstalk, verbose=True):
    det_xtalk = CrosstalkMatrix()
    for infile in file_list:
        if verbose:
            print "processing", infile
        ccd = MaskedCCD(infile, mask_files=mask_files)
        agg_amp, max_dn = aggressor(ccd)
        ratios = extractor(ccd, agg_amp)
        det_xtalk.set_row(agg_amp, ratios)
    return det_xtalk

if __name__ == '__main__':
    sys_xtfile = lambda amp : '/nfs/farm/g/lsst/u1/testData/eotestData/System/xtalk/data/xtalk_seg%02i.fits' % amp
    mask_files = ('CCD250_DEFECTS_mask.fits', )
    #
    # System crosstalk calculation
    #
    tstart = time.time()
    sys_xtalk = CrosstalkMatrix()
    for agg_amp in imutils.allAmps:
        ccd = MaskedCCD(sys_xtfile(agg_amp), mask_files=mask_files)
        ratios = system_crosstalk(ccd, agg_amp)
        sys_xtalk.set_row(agg_amp, ratios)
    print time.time() - tstart
    sys_xtalk.write('sys_xtalk.txt')
    #
    # Read it back in from the text file and plot.
    #
    foo = CrosstalkMatrix('sys_xtalk.txt')
    foo.plot_matrix('System crosstalk')

    #
    # Compute detector crosstalk from spot image datasets. (Use
    # system file as proxy.)
    #
    tstart = time.time()
    det_xtalk = CrosstalkMatrix()
    for agg_amp in imutils.allAmps:
        ccd = MaskedCCD(sys_xtfile(agg_amp), mask_files=mask_files)
        det_ratios = detector_crosstalk(ccd, agg_amp)
        det_xtalk.set_row(agg_amp, det_ratios)
    print time.time() - tstart

    tstart = time.time()
    det_xtalk_2 = CrosstalkMatrix()
    for agg_amp in imutils.allAmps:
        ccd = MaskedCCD(sys_xtfile(agg_amp), mask_files=mask_files)
        det_ratios = detector_crosstalk(ccd, agg_amp,
                                        signal_extractor=extract_mean_signal_2)
        det_xtalk_2.set_row(agg_amp, det_ratios)
    print time.time() - tstart

    sys_diff = sys_xtalk - det_xtalk
    print max(sys_diff.matrix.flat)
    print min(sys_diff.matrix.flat)

    sys_diff_2 = sys_xtalk - det_xtalk_2
    print max(sys_diff_2.matrix.flat)
    print min(sys_diff_2.matrix.flat)

    det_diff = det_xtalk - det_xtalk_2
    print max(det_diff.matrix.flat)
    print min(det_diff.matrix.flat)
