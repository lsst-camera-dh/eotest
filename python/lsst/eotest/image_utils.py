"""
@brief Module to perform standard operations on sensor images such
computing median images, unbiasing using the serial overscan region,
trimming, etc..
"""
import numpy as np
import numpy.random as random
import astropy.io.fits as fits
from fitsTools import fitsWriteto
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.pex.exceptions as pexExcept

class Metadata(object):
    def __init__(self, infile, hdu):
        self.header = None
        try:
            self.md = afwImage.readMetadata(infile, hdu)
        except:
            # This exception occurs when DM stack encounters a "." in
            # a FITS header keyword.
            self.header = fits.open(infile)[hdu-1].header
    def get(self, key):
        return self(key)
    def __call__(self, key):
        if self.header is None:
            return self.md.get(key)
        else:
            return self.header[key]

def allAmps(fits_file=None):
    all_amps = range(1, 17)
    if fits_file is None:
        return all_amps
    try:
        namps = fits.open(fits_file)[0].header['NAMPS']
        return range(1, namps+1)
    except KeyError:
        return all_amps

# Segment ID to HDU number in FITS dictionary
hdu_dict = dict( [ (1,'Segment10'), (2,'Segment11'), (3,'Segment12'),
                   (4,'Segment13'), (5,'Segment14'), (6,'Segment15'),
                   (7,'Segment16'), (8,'Segment17'), (9,'Segment07'),
                   (10,'Segment06'), (11,'Segment05'), (12,'Segment04'),
                   (13,'Segment03'), (14,'Segment02'), (15,'Segment01'),
                   (16,'Segment00') ] )

channelIds = dict([(i, hdu_dict[i][-2:]) for i in allAmps()])

mean = lambda x : afwMath.makeStatistics(x, afwMath.MEAN).getValue()
median = lambda x : afwMath.makeStatistics(x, afwMath.MEDIAN).getValue()
stdev = lambda x : afwMath.makeStatistics(x, afwMath.STDEV).getValue()

def dm_hdu(hdu):
    """ Compute DM HDU from the actual FITS file HDU."""
    return hdu + 1

def bias(im, overscan):
    """Compute the bias from the mean of the pixels in the serial
    overscan region."""
    return mean(im.Factory(im, overscan))

def bias_row(im, overscan, dxmin=5, dxmax=2):
    """Compute the mean signal for each row in the overscan region for
    a given amplifier on the CCD, skipping the first dxmin columns."""
    try:
        imarr = im.Factory(im, overscan).getArray()
    except AttributeError: # Dealing with a MaskedImage
        imarr = im.Factory(im, overscan).getImage().getArray()
    ny, nx = imarr.shape
    rows = np.arange(ny)
    values = np.array([np.mean(imarr[ii][dxmin:-dxmax]) for ii in rows])
    return(values)

def bias_func(im, overscan, fit_order=1, statistic=np.mean,
              nskip_cols=5, num_cols=15):
    """Compute the bias by fitting a polynomial (linear, by default)
    to the mean of each row of the selected overscan region.  This
    returns a numpy.poly1d object that returns the fitted bias as
    function of pixel row."""
    try:
        imarr = im.Factory(im, overscan).getArray()
    except AttributeError: # Dealing with a MaskedImage
        imarr = im.Factory(im, overscan).getImage().getArray()
    ny, nx = imarr.shape
    rows = np.arange(ny)
    if fit_order >= 0:
        values = np.array([statistic(imarr[j]) for j in rows])
        return np.poly1d(np.polyfit(rows, values, fit_order))
    else:
        # Use row-by-row bias level estimate, skipping initial columns
        # to avoid possible trailed charge.
        values = np.array([np.median(imarr[j][nskip_cols:nskip_cols+num_cols])
                           for j in rows])
        return lambda x: values[int(x)]

def bias_image(im, overscan, bias_method, fit_order=1, statistic=np.mean):
    biasim = afwImage.ImageF(im.getDimensions())
    imarr = biasim.getArray()
    ny, nx = imarr.shape
    if bias_method == 'bias_mean':
        my_bias = bias(im, overscan)
    elif bias_method == 'bias_row':
        my_bias = bias_row(im, overscan)
    elif bias_method == 'bias_func':
        poly = bias_func(im, overscan, fit_order, statistic)
        my_bias = poly(np.arange(ny))
    biasim = afwImage.ImageF(im.getDimensions())
    imarr = biasim.getArray()
    ny, nx = imarr.shape
    if isinstance(my_bias, float):
        my_bias = np.full(ny, my_bias)
    for row in range(ny):
        imarr[row] += my_bias[row]
    return biasim

def trim(im, imaging):
    "Trim the prescan and overscan regions."
    return im.Factory(im, imaging)

def unbias_and_trim(im, overscan, bias_method, fit_order=1, statistic=np.mean, 
                    imaging):
    """Subtract bias calculated from overscan region and optionally trim 
    prescan and overscan regions."""
    im -= bias_image(im, overscan, bias_method, fit_order, statistic)
    if imaging is not None:
        return trim(im, imaging)
    return im

def set_bitpix(hdu, bitpix):
    dtypes = {16 : np.int16, -32 : np.float32}
    for keyword in 'BSCALE BZERO'.split():
        if keyword in hdu.header.keys():
            del hdu.header[keyword]
    if bitpix > 0:
        my_round = np.round
    else:
        my_round = lambda x : x
    hdu.data = np.array(my_round(hdu.data), dtype=dtypes[bitpix])

def fits_median_file(files, outfile, bitpix=16, clobber=True):
    output = fits.open(files[0])
    for amp in allAmps(files[0]):
        try:
            del output[amp].header['BSCALE']
            del output[amp].header['BZERO']
        except KeyError:
            pass
        output[amp].data = fits_median(files, hdu=dm_hdu(amp)).getArray()
        if bitpix is not None:
            set_bitpix(output[amp], bitpix)
    fitsWriteto(output, outfile, clobber=clobber)

def fits_mean_file(files, outfile, bitpix=16, clobber=True):
    output = fits.open(files[0])
    all_amps = allAmps(files[0])
    for amp in all_amps:
        try:
            del output[amp].header['BSCALE']
            del output[amp].header['BZERO']
        except KeyError:
            pass
        output[amp].data = np.zeros(output[amp].data.shape)
    for infile in files:
        input = fits.open(infile)
        for amp in all_amps:
            output[amp].data += input[amp].data
    for amp in all_amps:
        output[amp].data /= len(files)
        if bitpix is not None:
            set_bitpix(output[amp], bitpix)
    fitsWriteto(output, outfile, clobber=clobber)

def fits_median(files, hdu=2, fix=True):
    """Compute the median image from a set of image FITS files."""
    ims = [afwImage.ImageF(f, hdu) for f in files]
    exptimes = [Metadata(f, 1).get('EXPTIME') for f in files]

    if min(exptimes) != max(exptimes):
        raise RuntimeError("Error: unequal exposure times")

    if fix:
        medians = np.array([median(im) for im in ims])
        med = sum(medians)/len(medians)
        errs = medians - med
        for im, err in zip(ims, errs):
            im -= err

    images = afwImage.vectorImageF()
    for image in ims:
        images.push_back(image)
    median_image = afwMath.statisticsStack(images, afwMath.MEDIAN)

    return median_image

def stack(ims, statistic=afwMath.MEDIAN):
    images = afwImage.vectorImageF()
    for image in ims:
        images.push_back(image)
    summary = afwMath.statisticsStack(images, statistic)
    return(summary)

def writeFits(images, outfile, template_file, bitpix=-32):
    output = fits.open(template_file)
    output[0].header['FILENAME'] = outfile
    for amp in images:
        output[amp].data = images[amp].getArray()
        set_bitpix(output[amp], bitpix)
    fitsWriteto(output, outfile, clobber=True, checksum=True)

def check_temperatures(files, tol, setpoint=None, warn_only=False):
    for infile in files:
        md = Metadata(infile, 1)  # Read PHDU keywords
        if setpoint is not None:
            ref_temp = setpoint
        else:
            ref_temp = md.get('TEMP_SET')
        ccd_temp = md.get('CCDTEMP')
        if np.abs(ccd_temp - ref_temp) > tol:
            what = "Measured operating temperature %(ccd_temp)s departs from expected temperature %(ref_temp)s by more than the %(tol)s tolerance for file %(infile)s" % locals()
            if warn_only:
                print what
            else:
                raise RuntimeError(what)

class SubRegionSampler(object):
    def __init__(self, dx, dy, nsamp, imaging):
        self.dx = dx
        self.dy = dy
        #
        # Generate sub-regions at random locations on the segment
        # imaging region.
        #
        self.xarr = random.randint(imaging.getWidth() - dx - 1, size=nsamp)
        self.yarr = random.randint(imaging.getHeight() - dy - 1, size=nsamp)
        self.imaging = imaging
    def bbox(self, x, y):
        return afwGeom.Box2I(afwGeom.Point2I(int(x), int(y)),
                             afwGeom.Extent2I(self.dx, self.dy))
    def subim(self, im, x, y):
        return im.Factory(im, self.bbox(x, y))

def bad_column(column_indices, threshold):
    """
    Count the sizes of contiguous sequences of masked pixels and
    return True if the length of any sequence exceeds the threshold
    number.
    """
    if len(column_indices) < threshold:
        # There are not enough masked pixels to mark this as a bad
        # column.
        return False
    # Fill an array with zeros, then fill with ones at mask locations.
    column = np.zeros(max(column_indices) + 1)
    column[(column_indices,)] = 1
    # Count pixels in contiguous masked sequences.
    masked_pixel_count = []
    last = 0
    for value in column:
        if value != 0 and last == 0:
            masked_pixel_count.append(1)
        elif value != 0 and last != 0:
            masked_pixel_count[-1] += 1
        last = value
    if len(masked_pixel_count) > 0 and max(masked_pixel_count) >= threshold:
        return True
    return False

def rebin_array(arr, binsize, use_mean=False):
    "See http://scipython.com/blog/binning-a-2d-array-in-numpy/"
    if binsize == 1:
        return arr
    shape = (arr.shape[0]//binsize, binsize, arr.shape[1]//binsize, binsize)
    if use_mean:
        result = arr.reshape(shape).mean(-1).mean(1)
    else:
        result = arr.reshape(shape).sum(-1).sum(1)
    return result

def rebin(image, binsize, use_mean=False):
    rebinned_array = rebin_array(image.getArray(), binsize, use_mean=use_mean)
    output_image = image.Factory(*rebinned_array.shape)
    out_array = output_image.getArray()
    out_array += rebinned_array
    return output_image

if __name__ == '__main__':
    import glob

    files = glob.glob('data/dark*.fits')
    im = fits_median(files)
