"""
@brief Module to perform standard operations on sensor images such
computing median images, unbiasing using the serial overscan region,
trimming, etc..
"""
from __future__ import print_function
from __future__ import absolute_import
from builtins import zip
from builtins import range
from builtins import object
import os
import warnings
import numpy as np
import numpy.random as random
from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning, AstropyUserWarning
from .fitsTools import fitsWriteto
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.pex.exceptions as pexExcept


class Metadata(object):
    def __init__(self, infile, hdu=0):
        self.header = None
        try:
            self.md = afwImage.readMetadata(infile, hdu)
        except:
            # This exception occurs when DM stack encounters a "." in
            # a FITS header keyword.
            self.header = fits.open(infile)[hdu].header

    def get(self, key):
        return self(key)

    def __call__(self, key):
        if self.header is None:
            return self.md.get(key)
        else:
            return self.header[key]


def allAmps(fits_file=None):
    all_amps = list(range(1, 17))
    if fits_file is None:
        return all_amps
    try:
        f = fits.open(fits_file)
        f.close()  # close first, in case exception is thrown and file is left open
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
    return hdu


def bias(im, overscan):
    """Compute the bias from the mean of the pixels in the serial
    overscan region."""
    return mean(im.Factory(im, overscan))


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


def bias_image(im, overscan, fit_order=1, statistic=np.mean):
    my_bias = bias_func(im, overscan, fit_order, statistic=statistic)
    biasim = afwImage.ImageF(im.getDimensions())
    imarr = biasim.getArray()
    ny, nx = imarr.shape
    for row in range(ny):
        imarr[row] += my_bias(row)
    return biasim


def trim(im, imaging):
    "Trim the prescan and overscan regions."
    return im.Factory(im, imaging)


def unbias_and_trim(im, overscan, imaging,
                    apply_trim=True, fit_order=1):
    """Subtract bias calculated from overscan region and optionally trim 
    prescan and overscan regions."""
    im -= bias_image(im, overscan, fit_order)
    if apply_trim:
        return trim(im, imaging)
    return im


def set_bitpix(hdu, bitpix):
    dtypes = {16: np.int16, -32: np.float32, 32: np.int32}
    for keyword in 'BSCALE BZERO'.split():
        if keyword in list(hdu.header.keys()):
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
        data  = fits_median(files, hdu=dm_hdu(amp)).getArray()
        output.append(fits.CompImageHDU(data=data,
                                        compresssion_type='RICE_1'))
        if bitpix is not None:
            set_bitpix(output[-1], bitpix)
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
                if bitpix is not None:
                    set_bitpix(output[amp], bitpix)
            fitsWriteto(output, outfile, overwrite=overwrite)


def fits_mean_file(files, outfile, overwrite=True, bitpix=32):
    output = fits.HDUList()
    output.append(fits.PrimaryHDU())
    all_amps = allAmps()
    for amp in all_amps:
        images = [afwImage.ImageF(item, amp) for item in files]
        mean_image = afwMath.statisticsStack(images, afwMath.MEAN)
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
                if bitpix is not None:
                    set_bitpix(output[amp], bitpix)
            for i in (-3, -2, -1):
                output.append(template[i])
            fitsWriteto(output, outfile, overwrite=overwrite)


def fits_median(files, hdu=1, fix=True):
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

    median_image = afwMath.statisticsStack(ims, afwMath.MEDIAN)

    return median_image


def writeFits(images, outfile, template_file, bitpix=-32):
    output = fits.HDUList()
    output.append(fits.PrimaryHDU())
    all_amps = allAmps()
    for amp in all_amps:
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
            for amp in all_amps:
                output[amp].header.update(template[amp].header)
                if bitpix is not None:
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
        ccd_temp = md.get('CCDTEMP')
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
    out_array += rebinned_array
    return output_image


if __name__ == '__main__':
    import glob

    files = glob.glob('data/dark*.fits')
    im = fits_median(files)
