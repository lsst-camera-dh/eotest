"""
@brief Class to handle masks for a full sensor.  Added mask files can be
or'd to existing masks, and mask bits can be set for use with an
afwMath.makeStatistics object.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import warnings
import astropy.io.fits as fits
from astropy.utils.exceptions import AstropyWarning, AstropyUserWarning
from lsst.eotest.fitsTools import fitsWriteto
from .AmplifierGeometry import makeAmplifierGeometry
from .ccd_bias_pca import CCD_bias_PCA
import lsst.daf.base as dafBase
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.ip.isr as ipIsr
import lsst.pex.exceptions as pexExcept
import lsst.eotest.image_utils as imutils

try:
    afwImage_Mask = afwImage.Mask
except AttributeError:
    afwImage_Mask = afwImage.MaskU


class MaskedCCDBiasImageException(RuntimeError):
    def __init__(self, *args):
        super(MaskedCCDBiasImageException, self).__init__(*args)


class MaskedCCD(dict):
    """
    This is the main abstraction for handling CCD data in the sensor
    acceptance test scripts.  The pixel data for each segment is
    represented by a MaskedImageF object and are accessed via the
    amplifier number.  Masks can be added and manipulated separately
    by various methods.
    """
    def __init__(self, imfile, mask_files=(), bias_frame=None,
                 interpolateFromMasks=False, linearity_correction=None,
                 dark_frame=None, all_amps=None):
        super(MaskedCCD, self).__init__()
        self.imfile = imfile
        self.md = imutils.Metadata(imfile)
        self.amp_geom = makeAmplifierGeometry(imfile)
        if all_amps is None:
            all_amps = imutils.allAmps(imfile)
        for amp in all_amps:
            image = afwImage.ImageF(imfile, imutils.dm_hdu(amp))
            mask = afwImage_Mask(image.getDimensions())
            self[amp] = afwImage.MaskedImageF(image, mask)
        self._added_mask_types = []
        for mask_file in mask_files:
            self.add_masks(mask_file)
        self.stat_ctrl = afwMath.StatisticsControl()
        if mask_files:
            self.setAllMasks()
        self.bias_frame = None
        self.ccd_pcas = None
        if bias_frame is not None:
            if isinstance(bias_frame, MaskedCCD):
                self.bias_frame = bias_frame
            elif isinstance(bias_frame, tuple) and len(bias_frame) == 2:
                # Assume this is a pair of filenames for the
                # CCD_bias_PCA model.
                self.ccd_pcas = CCD_bias_PCA.read_model(*bias_frame)
            else:
                self.bias_frame = MaskedCCD(bias_frame, all_amps=all_amps)
        if dark_frame is not None:
            self.dark_frame = MaskedCCD(dark_frame, bias_frame=bias_frame,
                                        all_amps=all_amps)
            self.dark_time_ratio \
                = self.md.get('DARKTIME')/self.dark_frame.md.get('DARKTIME')
        else:
            self.dark_frame = None
        self._interpolateFromMasks = interpolateFromMasks
        self._linearity_correction = linearity_correction

    def applyInterpolateFromMask(self, maskedImage, fwhm=0.001):
        try:
            for maskName in self._added_mask_types:
                try:
                    ipIsr.interpolateFromMask(maskedImage, fwhm=fwhm,
                                              maskName=maskName)
                except pexExcept.InvalidParameterError:
                    pass
        except TypeError:
            # interface change in ip_isr/isrFunctions.py v18.1.0:
            ipIsr.interpolateFromMask(maskedImage, fwhm=fwhm,
                                      maskNameList=self._added_mask_types)

    def mask_plane_dict(self):
        amp = list(self.keys())[0]
        return dict(list(self[amp].getMask().getMaskPlaneDict().items()))

    def add_masks(self, mask_file):
        """
        Add a masks from a mask file by or-ing with existing masks.
        """
        md = imutils.Metadata(mask_file)
        self._added_mask_types.append(md('MASKTYPE'))
        for amp in self:
            curr_mask = self[amp].getMask()
            curr_mask |= afwImage_Mask(mask_file, imutils.dm_hdu(amp))

    def setMask(self, mask_name=None, clear=False):
        """
        Enable a mask and return the afwMath.StatisticsControl object
        for use by afwMath.makeStatistics. If clear is False, then the
        new mask is or'd with the existing mask.  If clear is False
        and mask_name is None, then all mask bits are cleared.
        """
        if clear:                         # Unset all masks.
            self.stat_ctrl.setAndMask(0)
        if mask_name is not None:         # Set the desired mask.
            new_mask = (self.stat_ctrl.getAndMask()
                        | afwImage_Mask.getPlaneBitMask(mask_name))
            self.stat_ctrl.setAndMask(new_mask)
        return self.stat_ctrl

    def setAllMasks(self):
        "Enable all masks."
        mpd = self.mask_plane_dict()
        mask_bits = 2**len(mpd) - 1
        self.stat_ctrl.setAndMask(mask_bits)
        return self.stat_ctrl

    def write_bias_subtracted_MEF(self, outfile, gains=None, overwrite=True):
        """
        Write a bias-subtracted MEF file with the same format as
        the original raw FITS file.

        Parameters
        ----------
        outfile: str
            Output filename.
        gains: dict [None]
            Gains to apply to pixel data.  If None, then pixel values are
            written as ADUs.
        overwrite: bool [True]
            Flag to overwrite an existing output file.
        """
        hdulist = fits.HDUList()
        with fits.open(self.imfile) as template:
            hdulist.append(template[0])
            hdulist[0].header['ORIGFILE'] = hdulist[0].header['FILENAME']
            hdulist[0].header['FILENAME'] = outfile
            for amp in self:
                imarr = self.bias_subtracted_image(amp).getImage().getArray()
                if gains is not None:
                    imarr *= gains[amp]
                hdulist.append(fits.CompImageHDU(data=imarr,
                                                 header=template[amp].header))
            with warnings.catch_warnings():
                for warning in (UserWarning, AstropyWarning,
                                AstropyUserWarning):
                    warnings.filterwarnings('ignore', category=warning,
                                            append=True)
                fitsWriteto(hdulist, outfile, overwrite=True)


    def bias_image_using_overscan(self, amp, overscan=None, **kwargs):
        """
        Generate a bias image containing offset values calculated from
        bias(), bias_row(), bias_func() or bias_spline(). The default bias method
        is set to bias_row() in image_utils.py. Keyword arguments can be passed
        depending on which bias method is used.

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
        """
        if overscan is None:
            overscan = self.amp_geom.serial_overscan
        try:
            return imutils.bias_image(self._deep_copy(amp), overscan=overscan, **kwargs)
        except pexExcept.LSST_RUNTIME_EXCEPTION as eobj:
            raise MaskedCCDBiasImageException("DM stack error generating bias "
                                              + "image from overscan region:\n"
                                              + str(eobj))

    def _deep_copy(self, amp):
        return self[amp].Factory(self[amp], deep=True)

    def bias_image(self, amp, overscan=None, **kwargs):
        """
        Use separately stored metadata to determine file-specified
        overscan region. If bias_frame is not given, then calculate
        the bias image using bias(), bias_row(), bias_func() or bias_spline().
        The default bias method is set to bias_row() in image_utils.py.
        Keyword arguments can be passed depending on which bias method is used.

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
        """
        if self.bias_frame is not None:
            #
            # Use bias frame, if available, instead of overscan region
            #
            # Return a deep copy.
            my_image = self.bias_frame[amp].getImage()
            return my_image.Factory(my_image, deep=True)
        return self.bias_image_using_overscan(amp, overscan=overscan, **kwargs)

    def bias_subtracted_image(self, amp, overscan=None, **kwargs):
        """
        Subtract a bias image to correct for the offset. A bias correction is also
        applied if a base_frame is passed. The bias image with the offset values is
        generated using either of the bias(), bias_row(), bias_func() or bias_spline()
        methods from image_utils.py. The default bias method is set to bias_row().
        Keyword arguments can be passed depending on which bias method is used.

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
        """
        # Make a local copy to process and return.
        my_image = self._deep_copy(amp)
        if self.bias_frame is not None:
            # Make a deep copy of the bias frame.
            bias = self.bias_frame[amp].Factory(self.bias_frame[amp], deep=True)
            # Subtract x-independent component using overscan.
            bias -= \
                self.bias_frame.bias_image_using_overscan(amp,
                                                          overscan=overscan, **kwargs)
            # Subtract x-independent component of image for this amp
            # using overscan.
            my_image -= \
                self.bias_image_using_overscan(amp, overscan=overscan, **kwargs)
            # Subtract structured, x-dependent part.
            my_image -= bias
        elif self.ccd_pcas is not None:
            my_image.getImage().array[:] \
                -= self.ccd_pcas.pca_bias_correction(amp, my_image.getImage().array)
        else:
            my_image -= self.bias_image(amp, overscan, **kwargs)

        # Subtract dark frame, scaling the dark current by darktime ratio.
        if self.dark_frame is not None:
            dark_image = self.dark_frame.bias_subtracted_image(amp)
            dark_image *= self.dark_time_ratio
            my_image -= dark_image

        # Apply any linearity correction.
        if self._linearity_correction is not None:
            my_image.getImage().array[:] \
                = self._linearity_correction(amp, my_image.getImage().array)

        return my_image

    def unbiased_and_trimmed_image(self, amp, overscan=None,
                                   imaging=None, **kwargs):
        """
        Return an offset-corrected image where the offset values generated using
        either of the bias(), bias_row(), bias_func() or bias_spline() methods from
        image_utils.py. The default bias method is set to bias_row(). Keyword arguments
        can be passed depending on which bias method is used.

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
        """
        unbiased_image = self.bias_subtracted_image(amp, overscan, **kwargs)
        if imaging is None:
            imaging = self.amp_geom.imaging
        mi = imutils.trim(unbiased_image, imaging)
        if self._interpolateFromMasks:
            self.applyInterpolateFromMask(mi)
        return mi


def add_mask_files(mask_files, outfile, overwrite=True):
    amp_list = imutils.allAmps(mask_files[0])
    masks = dict([(amp, afwImage_Mask(mask_files[0], imutils.dm_hdu(amp)))
                  for amp in amp_list])
    for mask_file in mask_files[1:]:
        for amp in masks:
            masks[amp] |= afwImage_Mask(mask_file, imutils.dm_hdu(amp))
    output = fits.HDUList()
    output.append(fits.PrimaryHDU())
    output[0].header['MASKTYPE'] = 'SUMMED_MASKS'
    fitsWriteto(output, outfile, overwrite=overwrite)
    for amp in masks:
        md = dafBase.PropertySet()
        md.set('EXTNAME', 'SEGMENT%s' % imutils.channelIds[amp])
        masks[amp].writeFits(outfile, md, 'a')
    return masks


def compute_stats(image, sctrl, weights=None):
    flags = afwMath.MEAN | afwMath.STDEV
    if weights is None:
        stats = afwMath.makeStatistics(image, flags, sctrl)
    else:
        stats = afwMath.makeStatistics(image, weights, flags, sctrl)
    return stats.getValue(afwMath.MEAN), stats.getValue(afwMath.STDEV)


class MaskedCCDWrapper:
    """
    Wrapper class for MaskedCCD to reduce the memory footprint by
    having only one amp in memory at a time.
    """
    def __init__(self, imfile, **kwds):
        self.imfile = imfile
        self.kwds = kwds
        if 'all_amps' in kwds:
            self.all_amps = kwds['all_amps']
        self.ccd = None
        self.md = imutils.Metadata(imfile)

    def unbiased_and_trimmed_image(self, amp, **kwds):
        if self.ccd is None or not amp in self.ccd:
            self.kwds['all_amps'] = (amp,)
            self.ccd = MaskedCCD(self.imfile, **self.kwds)
        return self.ccd.unbiased_and_trimmed_image(amp, **kwds)

    def bias_subtracted_image(self, amp, **kwds):
        if self.ccd is None or not amp in self.ccd:
            self.kwds['all_amps'] = (amp,)
            self.ccd = MaskedCCD(self.imfile, **self.kwds)
        return self.ccd.bias_subtracted_image(amp, **kwds)

    def __getitem__(self, amp):
        if self.ccd is None or not amp in self.ccd:
            self.kwds['all_amps'] = (amp,)
            self.ccd = MaskedCCD(self.imfile, **self.kwds)
        return self.ccd[amp]

    def __getattr__(self, attr):
        if self.ccd is None:
            raise RuntimeError('__getattr__: self.ccd is None')
        return getattr(self.ccd, attr)


if __name__ == '__main__':
    image_file = 'bright_pix_test.fits'
    mask_files = ('bright_pix_mask.fits', 'CCD250_DEFECTS_mask.fits')

    ccd = MaskedCCD(image_file)
    for mask_file in mask_files:
        print("adding masks from", mask_file)
        ccd.add_masks(mask_file)
        print("mask plane dict:", ccd.mask_plane_dict())
        print()

    amp = imutils.allAmps()[0]

    sctrl = ccd.stat_ctrl
    print(sctrl.getAndMask(), compute_stats(ccd[amp], sctrl))

    sctrl = ccd.setMask('BAD')
    print(sctrl.getAndMask(), compute_stats(ccd[amp], sctrl))

    sctrl = ccd.setMask('CCD250_DEFECTS')
    print(sctrl.getAndMask(), compute_stats(ccd[amp], sctrl))

    sctrl = ccd.setMask(clear=True)
    print(sctrl.getAndMask(), compute_stats(ccd[amp], sctrl))

    sctrl = ccd.setAllMasks()
    print(sctrl.getAndMask(), compute_stats(ccd[amp], sctrl))
