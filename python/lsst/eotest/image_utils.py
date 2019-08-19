"""
@brief Module to perform standard operations on sensor images such
computing median images, unbiasing using the serial overscan region,
trimming, etc..
"""
from __future__ import print_function
from __future__ import absolute_import
import os
import warnings
import numpy as np
import numpy.random as random
from scipy import interpolate
from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning, AstropyUserWarning
from .fitsTools import fitsWriteto
import lsst.afw
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath

class Metadata(object):
    def __init__(self, infile, hdu=0):
        self.header = None
        try:
            self.md = afwImage.readMetadata(infile, dm_hdu(hdu))
        except:
            # This exception occurs when DM stack encounters a "." in
            # a FITS header keyword.
            self.header = fits.open(infile)[hdu].header

    def get(self, key):
        return self(key)

    def __call__(self, key):
        if self.header is None:
            return self.md.getScalar(key)
        else:
            return self.header[key]


def allAmps(fits_file=None):
    all_amps = list(range(1, 17))
    if fits_file is None:
        return all_amps
    try:
        with fits.open(fits_file) as f:
            namps = f[0].header['NAMPS']
        return list(range(1, namps+1))
    except KeyError:
        return all_amps


# Segment ID to HDU number in FITS dictionary
hdu_dict = dict([(1, 'Segment10'), (2, 'Segment11'), (3, 'Segment12'),
                 (4, 'Segment13'), (5, 'Segment14'), (6, 'Segment15'),
                 (7, 'Segment16'), (8, 'Segment17'), (9, 'Segment07'),
                 (10, 'Segment06'), (11, 'Segment05'), (12, 'Segment04'),
                 (13, 'Segment03'), (14, 'Segment02'), (15, 'Segment01'),
                 (16, 'Segment00')])

channelIds = dict([(i, hdu_dict[i][-2:]) for i in allAmps()])


def mean(x): return afwMath.makeStatistics(x, afwMath.MEAN).getValue()


def median(x): return afwMath.makeStatistics(x, afwMath.MEDIAN).getValue()


def stdev(x): return afwMath.makeStatistics(x, afwMath.STDEV).getValue()


def dm_hdu(hdu):
    """ Compute DM HDU from the actual FITS file HDU."""
    if lsst.afw.__version__.startswith('12.0'):
        return hdu + 1
    return hdu


def bias(im, overscan, **kwargs):
    """Compute the offset from the mean of the pixels in the serial
    overscan region.

    Args:
        im: A masked (lsst.afw.image.imageLib.MaskedImageF) or unmasked 
            (lsst.afw.image.imageLib.ImageF) afw image.
        overscan: A bounding box for the serial overscan region.

    Returns:
        A single float value for the mean of the overscan region.
    """
    return mean(im.Factory(im, overscan))

def bias_row(im, overscan, dxmin=5, dxmax=2, statistic=np.mean, **kwargs):
    """Compute the offset based on a statistic for each row in the serial 
    overscan region for columns dxmin through dxmax.

    Args:
        im: A masked (lsst.afw.image.imageLib.MaskedImageF) or unmasked 
            (lsst.afw.image.imageLib.ImageF) afw image.
        overscan: A bounding box for the serial overscan region.
        dxmin: The number of columns to skip at the beginning of the serial 
            overscan region.
        dxmax: The number of columns to skip at the end of the serial overscan region.
        statistic: The statistic to use to calculate the offset for each row.

    Returns:
        A numpy array with length equal to the number of rows in the serial overscan
        region.
    """
    try:
        imarr = im.Factory(im, overscan).getArray()
    except AttributeError: # Dealing with a MaskedImage
        imarr = im.Factory(im, overscan).getImage().getArray()
    ny, nx = imarr.shape
    rows = np.arange(ny)
    values = np.array([statistic(imarr[j][dxmin:-dxmax]) for j in rows])
    return lambda x: values[x]

def bias_func(im, overscan, dxmin=5, dxmax=2, statistic=np.mean, **kwargs):
    """Compute the offset by fitting a polynomial (order 1 by default)
    to the mean of each row of the serial overscan region.  This
    returns a numpy.poly1d object with the fitted bias as function of pixel row.
    Allows the option to explicitly set the fit order to apply to
    each row using additional **kwargs.

    Args:
        im: A masked (lsst.afw.image.imageLib.MaskedImageF) or unmasked
            (lsst.afw.image.imageLib.ImageF) afw image.
        overscan: A bounding box for the serial overscan region.
        dxmin: The number of columns to skip at the beginning of the serial
            overscan region.
        dxmax: The number of columns to skip at the end of the serial overscan region.
        statistic: The statistic to use to calculate the offset for each row.

    Keyword Arguments:
        fit_order: The order of the polynomial. The default is: 1.

    Returns:
        A np.poly1d object containing the coefficients for the polynomial fit.
    """
    try:
        imarr = im.Factory(im, overscan).getArray()
    except AttributeError: # Dealing with a MaskedImage
        imarr = im.Factory(im, overscan).getImage().getArray()
    ny, nx = imarr.shape
    rows = np.arange(ny)
    values = np.array([statistic(imarr[j][dxmin:-dxmax]) for j in rows])
    return np.poly1d(np.polyfit(rows, values, kwargs.get('fit_order', 1)))

def bias_spline(im, overscan, dxmin=5, dxmax=2, statistic=np.mean, **kwargs):
    """Compute the offset by fitting a spline to the mean of each row in the
    serial overscan region.

    Args:
        im: A masked (lsst.afw.image.imageLib.MaskedImageF) or unmasked
            (lsst.afw.image.imageLib.ImageF) afw image.
        overscan: A bounding box for the serial overscan region.
        dxmin: The number of columns to skip at the beginning of the serial
            overscan region.
        dxmax: The number of columns to skip at the end of the serial overscan region.
        statistic: The statistic to use to calculate the offset for each row.

    Keyword Arguments:
        k: The degree of the spline fit. The default is: 3.
        s: The amount of smoothing to be applied to the fit. The default is: 18000.
        t: The number of knots. If None, finds the number of knots to use
            for a given smoothing factor, s. The default is: None.

    Returns:
        A tuple (t,c,k) containing the vector of knots, the B-spline coefficients,
        and the degree of the spline.
    """

    try:
        imarr = im.Factory(im, overscan).getArray()
    except AttributeError: # Dealing with a MaskedImage
        imarr = im.Factory(im, overscan).getImage().getArray()
    ny, nx = imarr.shape
    rows = np.arange(ny)
    values = np.array([statistic(imarr[j][dxmin:-dxmax]) for j in rows])
    rms = 7 # Expected read noise per pixel
    weights = np.ones(ny) * (rms / np.sqrt(nx))
    return interpolate.splrep(rows, values, w=1/weights, k=kwargs.get('k', 3),
                              s=kwargs.get('s', 18000), t=kwargs.get('t', None))

def bias_image(im, overscan, dxmin=5, dxmax=2, statistic=np.mean, bias_method='row', **kwargs):
    """Generate a bias image containing the offset values calculated from
    bias(), bias_row(), bias_func() or bias_spline().

    Args:
        im: A masked (lsst.afw.image.imageLib.MaskedImageF) or unmasked
            (lsst.afw.image.imageLib.ImageF) afw image.
        overscan: A bounding box for the serial overscan region.
        dxmin: The number of columns to skip at the beginning of the serial
            overscan region.
        dxmax: The number of columns to skip at the end of the serial overscan region.
        statistic: The statistic to use to calculate the offset for each row.
        bias_method: Either 'mean', 'row', 'func' or 'spline'.

    Keyword Arguments:
        fit_order: The order of the polynomial. This only needs to be specified when
            using the 'func' method. The default is: 1.
        k: The degree of the spline fit. This only needs to be specified when using
            the 'spline' method. The default is: 3.
        s: The amount of smoothing to be applied to the fit. This only needs to be
            specified when using the 'spline' method. The default is: 18000.
        t: The number of knots. If None, finds the number of knots to use for a given
            smoothing factor, s. This only needs to be specified when using the 'spline'
            method. The default is: None.

    Returns:
        An image with size equal to the input image containing the offset level.
    """
    if bias_method not in ['mean', 'row', 'func', 'spline', 'none']:
        raise RuntimeError('Bias method must be either "none", "mean", "row", "func" or "spline".')  

    def dummy_none(im, overscan, dxmin, dxmax, **kwargs):
        return 0.0
    method = {'mean' : bias, 'row' : bias_row, 'func' : bias_func, 'spline' : bias_spline, 'none' : dummy_none}
    if bias_method not in ['mean', 'row', 'func', 'spline']:
        raise RuntimeError('Bias method must be either "mean", "row", "func" or "spline".')
    method = {'mean' : bias, 'row' : bias_row, 'func' : bias_func, 'spline' : bias_spline}
    my_bias = method[bias_method](im, overscan, dxmin=dxmin, dxmax=dxmax, **kwargs)
    biasim = afwImage.ImageF(im.getDimensions())
    imarr = biasim.getArray()
    ny, nx = imarr.shape
    if (bias_method == 'row') or (bias_method == 'func'):
        values = my_bias(np.arange(ny))
    elif bias_method == 'spline':
        values = interpolate.splev(np.arange(ny), my_bias)
    elif isinstance(my_bias, float):
        values = np.full(ny, my_bias)
    for row in range(ny):
        imarr[row] += values[row]
    biasim.setXY0(im.getX0(), im.getY0())
    return biasim

def trim(im, imaging):
    """Trim the prescan and overscan regions.

    Args:
        im: A masked (lsst.afw.image.imageLib.MaskedImageF) or unmasked
            (lsst.afw.image.imageLib.ImageF) afw image.
        imaging: A bounding box containing only the imaging section and
            excluding the prescan.
    Returns:
        An afw image.
    """

    return im.Factory(im, imaging)

def unbias_and_trim(im, overscan, imaging=None, dxmin=5, dxmax=2, bias_method='row',
                    bias_frame=None, **kwargs):
    """Subtract the offset calculated from the serial overscan region and optionally trim
    prescan and overscan regions. Includes the option to subtract the median of a stack of
    offset-subtracted bias frames to remove the bias level.

    Args:
        im: A masked (lsst.afw.image.imageLib.MaskedImageF) or unmasked
            (lsst.afw.image.imageLib.ImageF) afw image.
        overscan: A bounding box for the serial overscan region.
        imaging: A bounding box containing only the imaging section and
            excluding the prescan.
        dxmin: The number of columns to skip at the beginning of the serial
            overscan region.
        dxmax: The number of columns to skip at the end of the serial overscan region.
        bias_method: Either 'mean', 'row', 'func' or 'spline'.
        bias_frame: A single bias image containing a set of stacked oversan-corrected
            and trimmed bias frames.

    Keyword Arguments:
        fit_order: The order of the polynomial. This only needs to be specified when using
            the 'func' method. The default is: 1.
        k: The degree of the spline fit. This only needs to be specified when using the 'spline'
            method. The default is: 3.
        s: The amount of smoothing to be applied to the fit. This only needs to be specified when
            using the 'spline' method. The default is: 18000.
        t: The number of knots. If None, finds the number of knots to use for a given smoothing
            factor, s. This only needs to be specified when using the 'spline' method.
            The default is: None.

    Returns:
        An afw image.
    """

    im -= bias_image(im, overscan, dxmin=dxmin, dxmax=dxmax, bias_method=bias_method, **kwargs)
    if bias_frame:
        im -= bias_frame
    if imaging is not None:
        return trim(im, imaging)
    return im


def set_bitpix(hdu, bitpix):
    dtypes = {16: np.int16, -32: np.float32, 32: np.int32}
    for keyword in 'BSCALE BZERO'.split():
        if keyword in hdu.header:
            del hdu.header[keyword]
    if bitpix > 0:
        my_round = np.round
    else:
        def my_round(x): return x
    hdu.data = np.array(my_round(hdu.data), dtype=dtypes[bitpix])


def fits_median_file(files, outfile, bitpix=16, overwrite=True):
    output = fits.HDUList()
    output.append(fits.PrimaryHDU())
    all_amps = allAmps()
    for amp in all_amps:
        data = fits_median(files, hdu=dm_hdu(amp)).getArray()
        if bitpix < 0:
            output.append(fits.ImageHDU(data=data))
        else:
            output.append(fits.CompImageHDU(data=data,
                                            compresssion_type='RICE_1'))
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=UserWarning, append=True)
        warnings.filterwarnings('ignore', category=AstropyWarning, append=True)
        warnings.filterwarnings('ignore', category=AstropyUserWarning,
                                append=True)
        with fits.open(files[0]) as template:
            output[0].header.update(template[0].header)
            output[0].header['FILENAME'] = os.path.basename(outfile)
            for amp in all_amps:
                output[amp].header.update(template[amp].header)
                set_bitpix(output[amp], bitpix)
            fitsWriteto(output, outfile, overwrite=overwrite)


def fits_mean_file(files, outfile, overwrite=True, bitpix=32):
    output = fits.HDUList()
    output.append(fits.PrimaryHDU())
    all_amps = allAmps()
    for amp in all_amps:
        images = [afwImage.ImageF(item, dm_hdu(amp)) for item in files]
        if lsst.afw.__version__.startswith('12.0'):
            images = afwImage.vectorImageF(images)
        mean_image = afwMath.statisticsStack(images, afwMath.MEAN)
        if bitpix < 0:
            output.append(fits.ImageHDU(data=mean_image.getArray()))
        else:
            output.append(fits.CompImageHDU(data=mean_image.getArray(),
                                            compression_type='RICE_1'))
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=UserWarning, append=True)
        warnings.filterwarnings('ignore', category=AstropyWarning, append=True)
        warnings.filterwarnings('ignore', category=AstropyUserWarning,
                                append=True)
        with fits.open(files[0]) as template:
            output[0].header.update(template[0].header)
            output[0].header['FILENAME'] = os.path.basename(outfile)
            for amp in all_amps:
                output[amp].header.update(template[amp].header)
                set_bitpix(output[amp], bitpix)
            for i in (-3, -2, -1):
                output.append(template[i])
            fitsWriteto(output, outfile, overwrite=overwrite)


def fits_median(files, hdu=2, fix=True):
    """Compute the median image from a set of image FITS files."""
    ims = [afwImage.ImageF(f, hdu) for f in files]
    exptimes = [Metadata(f).get('EXPTIME') for f in files]

    if min(exptimes) != max(exptimes):
        raise RuntimeError("Error: unequal exposure times")

    if fix:
        medians = np.array([median(im) for im in ims])
        med = sum(medians)/len(medians)
        errs = medians - med
        for im, err in zip(ims, errs):
            im -= err

    if lsst.afw.__version__.startswith('12.0'):
        ims = afwImage.vectorImageF(ims)
    median_image = afwMath.statisticsStack(ims, afwMath.MEDIAN)

    return median_image

def stack(ims, statistic=afwMath.MEDIAN):
    """Stacks a list of images based on a statistic."""
    images = []
    for image in ims:
        images.append(image)
    if lsst.afw.__version__.startswith('12.0'):
        images = afwImage.vectorImageF(images)
    summary = afwMath.statisticsStack(images, statistic)
    return summary


def superbias(files, overscan, imaging=None, dxmin=5, dxmax=2, bias_method='row',
               hdu=2, statistic=afwMath.MEDIAN, **kwargs):
    """Generates a single stacked 'super' bias frame based on
    a statistic. Images must be either all masked or all unmasked."""
    ims = [afwImage.ImageF(f, hdu) for f in files]
    bias_frames = [unbias_and_trim(im, overscan, imaging, dxmin, dxmax, bias_method,
                                   **kwargs) for im in ims]
    return stack(bias_frames, statistic)

def superbias_file(files, overscan, outfile, imaging=None, dxmin=5, dxmax=2,
                    bias_method='row', bitpix=-32, clobber=True, **kwargs):
    images = {amp : superbias(files, overscan, imaging, dxmin, dxmax, bias_method,
                              hdu=dm_hdu(amp), **kwargs) for amp in allAmps(files[0])}
    writeFits(images, outfile, files[0], bitpix=bitpix)

def writeFits(images, outfile, template_file, bitpix=32):
    output = fits.HDUList()
    output.append(fits.PrimaryHDU())
    all_amps = allAmps()
    for amp in all_amps:
        if bitpix < 0:
            output.append(fits.ImageHDU(data=images[amp].getArray()))
        else:
            output.append(fits.CompImageHDU(data=images[amp].getArray(),
                                            compression_type='RICE_1'))
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=UserWarning, append=True)
        warnings.filterwarnings('ignore', category=AstropyWarning, append=True)
        warnings.filterwarnings('ignore', category=AstropyUserWarning,
                                append=True)

        with fits.open(template_file) as template:
            output[0].header.update(template[0].header)
            output[0].header['FILENAME'] = outfile
            metadata = images.get('METADATA', None)
            if metadata is not None:
                for key, val in metadata.items():
                    output[0].header[key] = val
            for amp in all_amps:
                output[amp].header.update(template[amp].header)
                set_bitpix(output[amp], bitpix)
            for i in (-3, -2, -1):
                output.append(template[i])
            fitsWriteto(output, outfile, overwrite=True, checksum=True)

def check_temperatures(files, tol, setpoint=None, warn_only=False):
    for infile in files:
        md = Metadata(infile)  # Read PHDU keywords
        if setpoint is not None:
            ref_temp = setpoint
        else:
            ref_temp = md.get('TEMP_SET')
        try:
            ccd_temp = md.get('CCDTEMP')
        except:
            print("Missing CCDTEMP keyword:", infile)
            return
        if np.abs(ccd_temp - ref_temp) > tol:
            what = "Measured operating temperature %(ccd_temp)s departs from expected temperature %(ref_temp)s by more than the %(tol)s tolerance for file %(infile)s" % locals(
            )
            if warn_only:
                print(what)
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
    out_array += rebinned_array.transpose()
    return output_image


if __name__ == '__main__':
    import glob

    files = glob.glob('data/dark*.fits')
    im = fits_median(files)
