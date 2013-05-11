import numpy as np
import lsst.afw.detection as afwDetect
import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath
import image_utils as imutils
from MaskedCCD import MaskedCCD
from BrightPixels import BrightPixels

def column_mean(raw_image, col, imaging=imutils.imaging):
    reg = afwGeom.Box2I(afwGeom.Point2I(col, imaging.getMinY()),
                        afwGeom.Extent2I(1, imaging.getHeight()))
    image = imutils.unbias_and_trim(raw_image, imaging=imaging)
    subim = image.Factory(image, reg)
    return afwMath.makeStatistics(subim, afwMath.MEAN).getValue()

def system_crosstalk(imfile, aggressor_amp, ethresh=10000,
                     mask_files=()):
    """
    Compute the system crosstalk.  ethresh is the threshold in DN for
    detecting the illuminated column in the aggressor amplifier.
    """
    #
    # Read in the system cross-talk images.
    #
    ccd = MaskedCCD(imfile, mask_files=mask_files)
    #
    # Find bright column of aggressor amplifier.
    #
    # Set exptime and gain to unity so that the e- threshold used by
    # BrightPixels converts directly to DN.
    #
    exptime = 1  
    gain = 1
    bp = BrightPixels(ccd[aggressor_amp], exptime=exptime, gain=gain,
                      ethresh=ethresh)
    pixels, columns = bp.find()

    if len(columns) > 1:
        raise RuntimeError("More than one aggressor column found.")

    agg_col = columns[0]
    agg_mean = column_mean(ccd[aggressor_amp], agg_col)
    
    ratios = dict([(amp, column_mean(ccd[amp], agg_col)/agg_mean)
                   for amp in imutils.allAmps])
    return ratios

def get_footprint(fp_set):
    footprints = [fp for fp in fp_set.getFootprints()]
    if len(footprints) > 1:
        raise RuntimeError("More than one spot image found in " +
                           "aggressor amplifier.")
    fp = footprints[0]
    peak_value = max([x.getPeakValue() for x in fp.getPeaks()])
    return fp, peak_value

def extract_mean_signal(masked_image, footprint):
    image = masked_image.getImage()
    maskarr = masked_image.getMask().getArray()
    image -= imutils.bias_image(image)
    imarr = image.getArray()
    signal = 0
    npix = 0
    for span in footprint.getSpans():
        y = span.getY()
        for x in range(span.getX0(), span.getX1()+1):
            npix += 1
            if maskarr[y][x] == 0:
                signal += imarr[y][x]
    return signal/float(npix)

def detector_crosstalk(imfile, aggressor_amp, ethresh=10000,
                       peak_frac=0.5, mask_files=()):
    """
    Compute detector crosstalk from a spot image in the aggressor
    amplifier.
    """
    ccd = MaskedCCD(imfile, mask_files=mask_files)
    image = imutils.unbias_and_trim(ccd[aggressor_amp])
    #
    # Extract footprint of spot image using nominal detection
    # threshold.
    #
    threshold = afwDetect.Threshold(ethresh)
    fp_set = afwDetect.FootprintSet(image, threshold)
    fp, peak_value = get_footprint(fp_set)
    #
    # Re-extract using threshold of peak_frac*peak_value to define the
    # "central region" of the spot image.
    #
    creg_threshold = afwDetect.Threshold(peak_frac*peak_value)
    creg_fp_set = afwDetect.FootprintSet(image, creg_threshold)
    footprint, peak_value = get_footprint(creg_fp_set)

    agg_mean = extract_mean_signal(ccd[aggressor_amp], footprint)
    ratios = dict([(amp, extract_mean_signal(ccd[amp], footprint)/agg_mean)
                   for amp in imutils.allAmps])
    return ratios

class CrosstalkMatrix(object):
    def __init__(self, filename=None, mode='w', namps=16):
        self.filename = filename
        self.mode = mode
        self.namps = namps
        if mode == 'r':
            self._read_matrix()
        else:
            self._set_matrix()
    def set_row(self, agg, ratios):
        self.matrix[agg-1] = np.array([ratios[amp] for amp in imutils.allAmps])
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
        output = open(outfile, 'w')
        output.write('#')
        for amp in range(1, self.namps+1, 1):
            output.write('victim%02i/aggressor  ' % amp)
        output.write('\n')
        for agg in range(self.namps):
            for victim in range(self.namps):
                output.write('%12.4e  ' % self.matrix[agg][victim])
            output.write('\n')
        output.close()
    def __sub__(self, other):
        result = CrosstalkMatrix()
        result.matrix = self.matrix - other.matrix
        return result
    def __add__(self, other):
        result = CrosstalkMatrix()
        result.matrix = self.matrix + other.matrix
        return result

if __name__ == '__main__':
    sys_xtfile = lambda amp : '/nfs/farm/g/lsst/u1/testData/eotestData/System/xtalk/data/xtalk_seg%02i.fits' % amp
    #
    # Loop over aggressor amps.
    #
    sys_xtalk = CrosstalkMatrix()
    det_xtalk = CrosstalkMatrix()
    for agg_amp in imutils.allAmps:
        ratios = system_crosstalk(sys_xtfile(agg_amp), agg_amp)
        sys_xtalk.set_row(agg_amp, ratios)

        det_ratios = detector_crosstalk(sys_xtfile(agg_amp), agg_amp)
        det_xtalk.set_row(agg_amp, det_ratios)

    print sys_xtalk.matrix
    print
    print det_xtalk.matrix
    print

    diff_xtalk = det_xtalk - sys_xtalk
    print diff_xtalk.matrix
