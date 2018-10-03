"""
@brief Crosstalk calculations for system and detector inputs.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import time
import numpy as np
import astropy.io.fits as fits
from lsst.eotest.fitsTools import fitsWriteto
import pylab
import pylab_plotter as plot
import lsst.afw.detection as afwDetect
import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath
import lsst.eotest.image_utils as imutils
from MaskedCCD import MaskedCCD
from BrightPixels import BrightPixels

def get_stats(image, stat_ctrl):
    flags = afwMath.MEDIAN | afwMath.STDEVCLIP
    stats = afwMath.makeStatistics(image, flags, stat_ctrl)
    return stats.getValue(afwMath.MEDIAN), stats.getValue(afwMath.STDEVCLIP)

def aggressor(ccd):
    """Guess the aggressor amp based on the maximum pixel value."""
    max_pix = lambda amp : max(ccd[amp].getImage().getArray().flat)
    candidate = ccd.keys()[0]
    max_pix_val = max_pix(candidate)
    for amp in ccd.keys()[1:]:
        val = max_pix(amp)
        if val > max_pix_val:
            candidate = amp
            max_pix_val = val
    return candidate, max_pix_val

def column_mean(ccd, amp, col, median_stack=None):
    imaging = ccd.amp_geom.imaging
    reg = afwGeom.Box2I(afwGeom.Point2I(col, imaging.getMinY()),
                        afwGeom.Extent2I(1, imaging.getHeight()))
    image = ccd.unbiased_and_trimmed_image(amp, median_stack=median_stack)
    subim = image.Factory(image, reg)
    flags = afwMath.MEAN | afwMath.STDEV | afwMath.NPOINT
    stats = afwMath.makeStatistics(subim, flags, ccd.stat_ctrl)
    mean = stats.getValue(afwMath.MEAN)
    npts = stats.getValue(afwMath.NPOINT)
    sigma = stats.getValue(afwMath.STDEV)/np.sqrt(npts)
    return np.array((mean, sigma))

def system_crosstalk(ccd, aggressor_amp, dnthresh=None, nsig=5, median_stack=None):
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
        image = ccd.unbiased_and_trimmed_image(aggressor_amp, median_stack=median_stack)
        median, stdev = get_stats(image, ccd.stat_ctrl)
        dnthresh = median + nsig*stdev
    bp = BrightPixels(ccd, aggressor_amp, exptime, gain, ethresh=dnthresh)
    pixels, columns = bp.find()

    if len(columns) > 1:
        raise RuntimeError("More than one aggressor column found.")

    agg_col = columns[0]
    agg_mean = column_mean(ccd, aggressor_amp, agg_coli, median_stack=median_stack)[0]

    ratios = dict([(amp, column_mean(ccd, amp, agg_col, median_stack=median_stack)/agg_mean)
                   for amp in ccd])
    return ratios

def get_footprint(fp_set, min_fp_size, threshold):
    footprints = [fp for fp in fp_set.getFootprints()
                  if fp.getNpix() >= min_fp_size]
    if len(footprints) > 1:
        message = "More than one spot image found in aggressor amplifier.\n"
        message += "      x     y     peak value  # pixels\n"
        for i, fp in enumerate(footprints):
            peak = [x for x in fp.getPeaks()][0]
            message += ('%2i  %4i  %4i     %6i       %4i\n' 
                        % (i, peak.getIx(), peak.getIy(), peak.getPeakValue(),
                           fp.getNpix()))
        message += "Threshold: %i\n" % threshold
        raise RuntimeError(message)
    fp = footprints[0]
    peak_value = max([x.getPeakValue() for x in fp.getPeaks()])
    return fp, peak_value

def extract_mean_signal_2(ccd, amp, footprint):
    masked_image = ccd.bias_subtracted_image(amp)
    stdev = afwMath.makeStatistics(masked_image, afwMath.STDEVCLIP,
                                   ccd.stat_ctrl).getValue()
    signal = 0
    npix = 0
    for span in footprint.getSpans():
        width = span.getX1() - span.getX0() + 1
        bbox = afwGeom.Box2I(afwGeom.Point2I(span.getX0(), span.getY()),
                             afwGeom.Extent2I(width, 1))
        subim = masked_image.Factory(masked_image, bbox)
        stats = afwMath.makeStatistics(subim, afwMath.SUM | afwMath.NPOINT,
                                       ccd.stat_ctrl)
        signal += stats.getValue(afwMath.SUM)
        npix += stats.getValue(afwMath.NPOINT)
    return np.array((signal/npix, stdev))

def extract_mean_signal(ccd, amp, footprint):
    maskarr = ccd[amp].getMask().getArray()
    image = ccd.bias_subtracted_image(amp)
    stdev = afwMath.makeStatistics(image, afwMath.STDEVCLIP,
                                   ccd.stat_ctrl).getValue()
    imarr = image.getImage().getArray()
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
                       signal_extractor=extract_mean_signal,
                       min_fp_size=50, median_stack=None):
    """
    Compute detector crosstalk from a spot image in the aggressor
    amplifier. dnthresh is the threshold in DN for detecting the
    illuminated column in the aggressor amplifier; if set to None,
    then nsig*clipped_stdev above median is used for the threshold.
    """
    image = ccd.unbiased_and_trimmed_image(aggressor_amp, median_stack=median_stack)
    #
    # Extract footprint of spot image using nominal detection
    # threshold.
    #
    if dnthresh is None:
        median, stdev = get_stats(image, ccd.stat_ctrl)
#        dnthresh = median + nsig*stdev
        dnthresh = (np.max(ccd[aggressor_amp].getImage().getArray())
                    + median)/2.
#    print "dnthresh =", dnthresh
    threshold = afwDetect.Threshold(dnthresh)
    fp_set = afwDetect.FootprintSet(image, threshold)
    try:
        footprint, peak_value = get_footprint(fp_set, min_fp_size, dnthresh)
    except IndexError:
        raise RuntimeError('index error in get_footprint')

    agg_mean = signal_extractor(ccd, aggressor_amp, footprint)[0]
    ratios = dict([(amp, signal_extractor(ccd, amp, footprint)
                    /agg_mean) for amp in ccd])
#    for amp in ratios:
#        if ratios[amp][0] > 0.1:
#            ratios[amp] = (0, 0)
    return ratios

class CrosstalkMatrix(object):
    def __init__(self, filename=None, namps=16):
        self.filename = filename
        self.namps = namps
        self._set_matrix()
        if self.filename is not None:
            self._read_matrix()
    def set_row(self, agg, ratios):
        self.matrix[agg-1] = np.array([ratios[amp][0] for amp
                                       in imutils.allAmps()])
    def _set_matrix(self):
        self.matrix = np.zeros((self.namps, self.namps), dtype=np.float)
    def _read_matrix(self):
        if self.filename[-5:] == '.fits':
            self._read_fits_matrix()
        else:
            self._read_text_matrix()
    def _read_fits_matrix(self):
        self.matrix = fits.open(self.filename)[0].data
    def _read_text_matrix(self):
        input = open(self.filename, 'r')
        amp = 0
        for line in input:
            if line[0] == '#':
                continue
            self.matrix[amp] = np.array([float(x) for x
                                         in line.strip().split()])
            amp += 1
    def write_fits(self, outfile=None, clobber=True):
        if outfile is None:
            outfile = self.filename
        else:
            self.filename = outfile
        output = fits.HDUList()
        output.append(fits.PrimaryHDU(data=self.matrix))
        fitsWriteto(output, outfile, clobber=clobber)
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
                    scale_factor=1e2, fontsize=10, figsize=(12, 6),
                    cmap=None):
        pylab.ion()
        cmin, cmax = cmap_range
        my_matrix = np.copy(self.matrix)*scale_factor
        for i in range(self.namps):
            my_matrix[i][i] = 0
        thresh = 10**(-precision)
        pylab.ion()
        if cmap is None:
            cdict = dict(red=((0, cmin, cmin), (1, cmax, cmax)),
                         green=((0, cmin, cmin), (1, cmax, cmax)),
                         blue=((0, cmin, cmin), (1, cmax, cmax)))
            cmap = pylab.matplotlib.colors.LinearSegmentedColormap('my_grey',
                                                                   cdict, 256)
        fig = pylab.figure(figsize=figsize)
        axes = fig.add_subplot(111)
        image = pylab.imshow(np.abs(my_matrix), interpolation='nearest',
                             aspect='auto', cmap=cmap)
        pylab.xlabel('victim')
        pylab.ylabel('aggressor')
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
        pylab.colorbar(image)
    def plot(self, cmap=pylab.cm.hot, title=''):
        my_matrix = np.copy(self.matrix)
        for i in range(self.namps):
            my_matrix[i][i] = 0
        win = plot.Window()
        fig, axes = win.fig, win.axes[-1]
        foo = axes.imshow(my_matrix, interpolation='nearest', cmap=cmap)
        pylab.xlabel('victim')
        pylab.ylabel('aggressor')
        axes.set_xticks(range(self.namps))
        axes.set_yticks(range(self.namps))
        axes.set_xticklabels(['%i' % i for i in range(1, self.namps+1)])
        axes.set_yticklabels(['%i' % i for i in range(1, self.namps+1)])
        cbar = fig.colorbar(foo)
        axes.set_title(title)
    def __sub__(self, other):
        result = CrosstalkMatrix()
        result.matrix = self.matrix - other.matrix
        return result
    def __add__(self, other):
        result = CrosstalkMatrix()
        result.matrix = self.matrix + other.matrix
        return result

def make_crosstalk_matrix(file_list, mask_files=(),
                          extractor=detector_crosstalk, verbose=True, median_stack=None):
    det_xtalk = CrosstalkMatrix()
    try:
        os.path.isfile(file_list)
        # A single file, so we assume that we have a multi-aggressor spot frame.
        ccd = MaskedCCD(file_list, mask_files=mask_files)
        for agg_amp in ccd:
            if verbose:
                print "processing aggressor amp", agg_amp
            try:
                det_ratios = extractor(ccd, agg_amp, median_stack=median_stack)
                det_xtalk.set_row(agg_amp, det_ratios)
            except RuntimeError, message:
                print "Error extracting victim/aggressor ratios."
                print message
                print "Skipping."
    except TypeError:
        # Presumably, we have a 16 amplifier dataset.
        for infile in file_list:
            if verbose:
                print "processing", infile
            ccd = MaskedCCD(infile, mask_files=mask_files)
            agg_amp, max_dn = aggressor(ccd)
            try:
                ratios = extractor(ccd, agg_amp, median_stack=median_stack)
                det_xtalk.set_row(agg_amp, ratios)
            except RuntimeError, message:
                print "Error extracting victim/aggressor ratios:"
                print message
                print "Skipping."
    return det_xtalk

if __name__ == '__main__':
    sys_xtfile = lambda amp : '/u/gl/jchiang/ki18/LSST/SensorTests/eotest/0.0.0.6/work/system/xtalk/debug/xtalk_%02i_debug.fits' % amp
#    mask_files = ('CCD250_DEFECTS_mask.fits', )
    mask_files = ()
    #
    # System crosstalk calculation
    #
    tstart = time.time()
    sys_xtalk = CrosstalkMatrix()
    for agg_amp in imutils.allAmps():
        ccd = MaskedCCD(sys_xtfile(agg_amp), mask_files=mask_files)
        ratios = system_crosstalk(ccd, agg_amp)
        sys_xtalk.set_row(agg_amp, ratios)
    print time.time() - tstart
    sys_xtalk.write('sys_xtalk.txt')
    sys_xtalk.write_fits('sys_xtalk.fits')
    #
    # Read it back in from the fits file and plot.
    #
#    foo = CrosstalkMatrix('sys_xtalk.txt')
    foo = CrosstalkMatrix('sys_xtalk.fits')
    foo.plot_matrix('System crosstalk')

    #
    # Compute detector crosstalk from spot image datasets. (Use
    # system file as proxy.)
    #
    tstart = time.time()
    det_xtalk = CrosstalkMatrix()
    for agg_amp in imutils.allAmps():
        ccd = MaskedCCD(sys_xtfile(agg_amp), mask_files=mask_files)
        det_ratios = detector_crosstalk(ccd, agg_amp)
        det_xtalk.set_row(agg_amp, det_ratios)
    print time.time() - tstart

    tstart = time.time()
    det_xtalk_2 = CrosstalkMatrix()
    for agg_amp in imutils.allAmps():
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
