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
from scipy.ndimage.filters import gaussian_filter
from scipy import optimize
from MaskedCCD import MaskedCCD
from BrightPixels import BrightPixels

import glob # just used for testing

def get_stats(image, stat_ctrl):
    flags = afwMath.MEDIAN | afwMath.STDEVCLIP
    stats = afwMath.makeStatistics(image, flags, stat_ctrl)
    return stats.getValue(afwMath.MEDIAN), stats.getValue(afwMath.STDEVCLIP)

def get_spot_array(ccd, amp, ay, ax, r=50):
    """Apply circular mask to capture a spot on the amplifier."""
    imarr = ccd[amp].getImage().getArray()
    ny, nx = imarr.shape
    
    y, x = np.ogrid[-ay:ny-ay, -ax:nx-ax]
    mask = x*x + y*y >= r*r
    
    spot_arr = np.ma.MaskedArray(imarr, mask)

    return spot_arr

def aggressor(ccd, pix_threshold=10000, sigma=50):
    """Locate aggressor spots on the CCD amplifiers."""

    candidate_list = []
    for amp in ccd.keys():
        blurred = gaussian_filter(ccd[amp].getImage().getArray(), sigma)
        y, x = np.unravel_index(blurred.argmax(), blurred.shape)
        candidate_spot = get_spot_array(ccd, amp, y, x, r=50)

        if np.mean(candidate_spot) > pix_threshold:
            candidate_list.append(amp)

    return candidate_list

def victim_model(var_array, aggressor):
    """Make a model of a cross talk victim postage stamp."""
    
    # Get coefficients
    xtalk_signal = var_array[0]
    bias = var_array[1]
    tilty = var_array[2]
    tiltx = var_array[3]

    Y, X = np.mgrid[:aggressor.shape[0], :aggressor.shape[1]]
    model = xtalk_signal*aggressor + tilty*Y + tiltx*X + bias

    return model

def crosstalk_model_fit(ccd, aggressor_amp, amp, footprint, num_iter=10):
    """Perform an iterative least-squares model fit for a cross talk victim."""

    bbox = footprint.getBBox()
    grow = (200-min(bbox.getWidth(), bbox.getHeight()))/2
    bbox.grow(grow)

    aggressor_image = ccd.bias_subtracted_image(aggressor_amp)
    sub_image = aggressor_image.getImage().Factory(aggressor_image.getImage(), bbox)
    aggressor_stamp = sub_image.getArray()

    victim_image = ccd.bias_subtracted_image(amp)
    sub_image = victim_image.getImage().Factory(victim_image.getImage(), bbox)
    victim_stamp = sub_image.getArray()

    crosstalk_results = np.asarray([[0,0,0,0]])
    victim_array = np.ma.masked_invalid(victim_stamp)
    mask = np.ma.getmask(victim_array)

    for i in range(num_iter):

        ## Mask outliers using residual
        model = np.ma.masked_where(mask, victim_model(crosstalk_results[0], 
                                                      aggressor_stamp))
        residual = victim_array - model
        victim_array = np.ma.masked_where(np.abs(residual-residual.mean()) \
                                          > 3.0*residual.std(), victim_stamp)
        mask = np.ma.getmask(victim_array)
                                        
        ## Construct basis arrays
        background = np.ma.masked_where(mask, np.ones(aggressor_stamp.shape))
        Y, X = np.mgrid[:aggressor_stamp.shape[0], :aggressor_stamp.shape[1]]
        Y = np.ma.masked_where(mask, Y)
        X = np.ma.masked_where(mask, X)
        aggressor_array = np.ma.masked_where(mask, aggressor_stamp)

        ## Perform least-square minimization                                 
        b = victim_array.compressed()/victim_array.std()
        A = np.vstack([aggressor_array.compressed(),
                       background.compressed(),
                       Y.compressed(), 
                       X.compressed()]).T/victim_array.std()

        crosstalk_results = np.linalg.lstsq(A, b)
        covar = np.matrix(np.dot(A.T, A)).I

    return np.append(crosstalk_results[0], covar[0,0])
        

def get_footprint(fp_set, min_fp_size, threshold):
    """Locate and return a cross talk spot footprint."""
    
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
        
def detector_crosstalk(ccd, aggressor_amp, dnthresh=None, nsig=5,
                       min_fp_size=50):
    """
    Compute detector crosstalk from a spot image in the aggressor
    amplifier. dnthresh is the threshold in DN for detecting the
    illuminated column in the aggressor amplifier; if set to None,
    then nsig*clipped_stdev above median is used for the threshold.
    """
    image = ccd.unbiased_and_trimmed_image(aggressor_amp)

    if dnthresh is None:
        median, stdev = get_stats(image, ccd.stat_ctrl)
        dnthresh = (np.max(ccd[aggressor_amp].getImage().getArray()) + median)/2.

    threshold = afwDetect.Threshold(dnthresh)
    fp_set = afwDetect.FootprintSet(image, threshold)
    try:
        footprint, peak_value = get_footprint(fp_set, min_fp_size, dnthresh)
    except IndexError:
        raise RuntimeError('index error in get_footprint')

    crosstalk_results = dict([(amp, crosstalk_model_fit(ccd, aggressor_amp, amp, footprint))
                             for amp in ccd])
    
    return crosstalk_results

class CrosstalkMatrix():
    
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

def make_crosstalk_matrix(file_list, mask_files=(), extractor=detector_crosstalk, verbose=True):
    
    det_xtalk = CrosstalkMatrix()
    for infile in file_list:
        if verbose:
            print "processing", infile
        ccd = MaskedCCD(infile, mask_files=mask_files)
        agg_amp_list = aggressor(ccd)
        for agg_amp in agg_amp_list:
            try:
                ratios = extractor(ccd, agg_amp)
                det_xtalk.set_row(agg_amp, ratios)
            except RuntimeError, message:
                print "Error extracting victim/aggressor ratios:"
                print message
                print "Skipping."           
    return det_xtalk

if __name__ == '__main__':

    test_files = glob.glob('/nfs/slac/g/ki/ki19/lsst/snyder18/LSST/Data/test/S11_00*_calibrated.fits')
    print test_files
    test_xtalk_result = make_crosstalk_matrix(test_files)
