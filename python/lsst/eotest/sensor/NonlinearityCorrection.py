"""
Code to apply non-linearity correction.
"""
import copy

import numpy as np

import scipy.optimize
from scipy.interpolate import UnivariateSpline

import astropy.io.fits as fits

from lsst.eotest.fitsTools import fitsTableFactory, fitsWriteto


def lin_func(pars, xvals):
    """Return a line whose slope is pars[0]"""
    return pars[0]*xvals

def chi2_model(pars, xvals, yvals):
    """Return the chi2 w.r.t. the model"""
    return (yvals - lin_func(pars, xvals))/np.sqrt(yvals)

def make_profile_hist(xbin_edges, xdata, ydata, **kwargs):
    """Build a profile historgram

    Parameters
    ----------
    xbin_edges : `array`
        The bin edges
    xdata : `array`
        The x-axis data
    ydata : `array`
        The y-axis data

    Keywords
    --------
    yerrs :  `array`
        The errors on the y-axis points

    stderr : `bool`
        Set error bars to standard error instead of RMS

    Returns
    -------
    x_vals : `array`
        The x-bin centers
    y_vals : `array`
        The y-bin values
    y_errs : `array`
        The y-bin errors
    """
    yerrs = kwargs.get('yerrs', None)
    stderr = kwargs.get('stderr', False)

    nx = len(xbin_edges) - 1
    x_vals = (xbin_edges[0:-1] + xbin_edges[1:])/2.
    y_vals = np.ndarray((nx))
    y_errs = np.ndarray((nx))

    if yerrs is None:
        weights = np.ones(ydata.shape)
    else:
        weights = 1./(yerrs*yerrs)

    y_w = ydata*weights

    for i, (xmin, xmax) in enumerate(zip(xbin_edges[0:-1], xbin_edges[1:])):
        mask = (xdata >= xmin) * (xdata < xmax)
        if mask.sum() < 2:
            y_vals[i] = 0.
            y_errs[i] = -1.
            continue
        y_vals[i] = y_w[mask].sum() / weights[mask].sum()
        y_errs[i] = ydata[mask].std()
        if stderr:
            y_errs[i] /= np.sqrt(mask.sum())

    return x_vals, y_vals, y_errs


class NonlinearityCorrection:
    """Apply a non-linearity correction

    The point of this calls is to serve as a callable object that will
    linearize bias-subtracted data

    corrected_adu = nlc(amp, uncorrected_adu)

    This is implemented as a spline interpolation for each of the 16 amplifiers on a CCD
    """
    def __init__(self, prof_x, prof_y, prof_yerr, **kwargs):
        """C'tor

        Parameters
        ----------
        prof_x : `array`
            Array of 16 x nbins values for the x-axis of the correction function
        prof_y : `array`
            Array of 16 x nbins values for the y-axis of the correction function
        prof_yerr : `array`
            Array of 16 x nbins values for the y-axis of the correction function

        Keywords
        --------
        Passed to UnivariateSpline c'tor

        """
        self._prof_x = prof_x
        self._prof_y = prof_y
        self._prof_yerr = prof_yerr
        self._nxbins = self._prof_x.shape[1]

        kwcopy = kwargs.copy()
        kwcopy.setdefault('s', 1e-6)
        kwcopy.setdefault('ext', 3)

        self._spline_dict = {}
        for iamp in range(16):
            idx_sort = np.argsort(self._prof_x[iamp])
            profile_x = self._prof_x[iamp][idx_sort]
            profile_y = self._prof_y[iamp][idx_sort]
            if self._prof_yerr is not None:
                profile_yerr = self._prof_yerr[iamp][idx_sort]
                mask = profile_yerr >= 0.
            else:
                mask = np.ones(profile_x.shape)
            try:
                self._spline_dict[iamp] = UnivariateSpline(profile_x[mask],
                                                           profile_y[mask],
                                                           **kwcopy)
            except Exception:
                self._spline_dict[iamp] = lambda x : x

    def __getitem__(self, amp):
        """Get the function that corrects a particular amp"""
        return self._spline_dict[amp]

    def __call__(self, amp, adu):
        """Apply the non-linearity correction to a particular amp"""
        return adu*self._spline_dict[amp-1](adu)


    def write_to_fits(self, fits_file):
        """Write this object to a FITS file"""
        output = fits.HDUList()
        output.append(fits.PrimaryHDU())

        col_prof_x = fits.Column(name='prof_x', format='%iE' % self._nxbins,
                                 unit='ADU', array=self._prof_x)
        col_prof_y = fits.Column(name='prof_y_corr', format='%iE' % self._nxbins,
                                 unit='ADU', array=self._prof_y)
        col_prof_yerr = fits.Column(name='prof_yerr', format='%iE' % self._nxbins,
                                    unit='ADU', array=self._prof_yerr)

        fits_cols = [col_prof_x, col_prof_y, col_prof_yerr]
        hdu = fitsTableFactory(fits_cols)
        hdu.name = 'nonlin'
        output.append(hdu)

        fitsWriteto(output, fits_file, overwrite=True)


    def save_plots(self, plotfile, **kwargs):
        """Save plots showing the nonlinearity correction"""
        import matplotlib.pyplot as plt
        ymin = kwargs.get('ymin', None)
        ymax = kwargs.get('ymax', None)

        figsize = kwargs.get('figsize', (15, 10))

        fig, axs = plt.subplots(nrows=4, ncols=4, figsize=figsize)
        fig.suptitle("Nonlinearity")

        xlabel = r'Mean [ADU]'
        ylabel = r'Frac Resid [$(q - g\mu)/g\mu$]'
        for i_row in range(4):
            ax_row = axs[i_row, 0]
            ax_row.set_ylabel(ylabel)

        for i_col in range(4):
            ax_col = axs[3, i_col]
            ax_col.set_xlabel(xlabel)

        iamp = 0
        for i_row in range(4):
            for i_col in range(4):
                axes = axs[i_row, i_col]
                if ymin is not None or ymax is not None:
                    axes.set_ylim(ymin, ymax)
                mask = self._prof_yerr[iamp] >= 0.
                x_masked = self._prof_x[iamp][mask]
                xline = np.linspace(1., x_masked.max(), 1001)
                model = self._spline_dict[iamp](xline)
                axes.errorbar(x_masked, self._prof_y[iamp][mask],
                              yerr=self._prof_yerr[iamp][mask], fmt='.')
                axes.plot(xline, model, 'r-')
                iamp += 1
        if plotfile is None:
            fig.show()
        else:
            fig.savefig(plotfile)


    @classmethod
    def create_from_table(cls, table, **kwargs):
        """Create a NonlinearityCorrection object from a fits file

        Parameters
        ----------
        table : `Table`
            The table data used to build the nonlinearity correction

        kwargs : passed to UnivariateSpline Constructor

        Returns
        -------
        nl : `NonlinearityCorrection`
            The requested object
        """
        prof_x = table.data['prof_x']
        prof_y = table.data['prof_y_corr']
        prof_yerr = table.data['prof_yerr']
        return cls(prof_x, prof_y, prof_yerr, **kwargs)

    @classmethod
    def create_from_fits_file(cls, fits_file, hdu_name='nonlin', **kwargs):
        """Create a NonlinearityCorrection object from a fits file

        Parameters
        ----------
        fits_file : `str`
            The file with the data used to build the nonlinearity correction

        hdu_name : `str`
            The name of the HDU with the nonlinearity correction data

        kwargs : passed to UnivariateSpline Constructor

        Returns
        -------
        nl : `NonlinearityCorrection`
            The requested object
        """
        hdulist = fits.open(fits_file)
        table = hdulist[hdu_name]
        nl = cls.create_from_table(table, **kwargs)
        hdulist.close()
        return nl


    @staticmethod
    def _correct_unit_point(profile_x, profile_y, profile_yerr, unit_point):
        """Force the spline to go through one at a particular x-value

        Parameters
        ----------
        profile_x : `array`
            The x-bin centers
        profile_y : `array`
            The b-bin values
        profile_yerr : `array`
            The y-bin errors
        unit_point : `float`
            The x-value where the spline should go through one

        Returns
        -------
        y_vals_corr
            The adjusted y-values
        y_errs_corr
            The adjusted y-errors
        """
        uni_spline = UnivariateSpline(profile_x, profile_y)
        y_value_0 = uni_spline(unit_point)

        y_vals_corr = profile_y/y_value_0
        y_errs_corr = profile_yerr
        return y_vals_corr, y_errs_corr

    @classmethod
    def create_from_det_response(cls, detresp, gains, **kwargs):
        """Create a NonlinearityCorrection object DetectorResponse FITS file

        Note that the DetectorResponse files typically store the signal in electrons,
        but we want a correction that works on ADU, so we have to remove the gains.

        Parameters
        ----------
        detresp : `DetectorResponse`
            An object with the detector response calculated from flat-pair files

        gains : `array` or `None`
            Array with amplifier by amplifer gains


        Keywords
        --------
        fit_range : `tuple`
            The range over which to define the non-linearity, defaults to (0., 9e4)

        nprofile_bins : `int` or `None`
             The number of bins to use in the profile, defaults to 10
             If `None` then this will use all of the data point rather that making
             a profile histogram

        null_point : `float` or `None`
             X-value at which the correction should vanish, defaults to 0.
             If `None` then this will simply use the pivot point of the fit to the data

        remaining kwargs are passed to the class c'tor

        Returns
        -------
        nl : `NonlinearityCorrection`
            The requested object
        """
        kwcopy = kwargs.copy()
        fit_range = kwcopy.pop('fit_range', (0., 9e4))
        unit_point = kwcopy.pop('unit_point', None)
        nprofile_bins = kwcopy.pop('nprofile_bins', 100)
        nprofile_bins = min(len(detresp.flux)//3, nprofile_bins)

        if nprofile_bins is not None:
            xbins = np.linspace(fit_range[0], fit_range[1], nprofile_bins+1)
        else:
            xbins = None
            nprofile_bins = len(detresp.flux)

        prof_x = np.ndarray((16, nprofile_bins))
        prof_y = np.ndarray((16, nprofile_bins))
        prof_yerr = np.ndarray((16, nprofile_bins))

        for idx, amp in enumerate(detresp.Ne):
            # For each amp, fit a linear model to the signal vs
            # incident flux.
            xdata = copy.copy(detresp.flux)
            ydata = copy.copy(detresp.Ne[amp])
            # The nominal fit_range applies to e-/pixel.  We also want
            # to avoid fitting data past the saturation peak.
            ypeak_index = np.argmax(ydata)
            x_at_ypeak = xdata[ypeak_index]
            mask = np.where((fit_range[0] < ydata) & (fit_range[1] > ydata)
                            & (xdata < x_at_ypeak))
            # Convert to ADU
            ydata /= gains[amp]

            xdata_fit = xdata[mask]
            ydata_fit = ydata[mask]
            mean_slope = (ydata_fit/xdata_fit).mean()
            pars = (mean_slope,)
            results = scipy.optimize.leastsq(chi2_model, pars,
                                             full_output=1,
                                             args=(xdata_fit, ydata_fit))
            model_yvals = lin_func(results[0], xdata)

            # Compute the ratio of the linear model to the measured
            # signal in ADU.  This ratio would be the correction
            # factor assuming the y value computed from the linear fit
            # is the desired signal.  The correction factor would then
            # be a function of the measured signal that's in ydata.
            ratio = model_yvals/ydata

            # Don't try to fit a spline past the saturation peak.
            ydata = ydata[:ypeak_index]
            ratio = ratio[:ypeak_index]

            if xbins is not None:
                prof_x[idx], prof_y[idx], prof_yerr[idx] \
                    = make_profile_hist(xbins, ydata, ratio, stderr=True)
            else:
                prof_x[idx], prof_y[idx], prof_yerr[idx] = ydata, ratio, None

            if unit_point is not None:
                prof_y[idx], prof_yerr[idx] \
                    = cls._correct_unit_point(prof_x[idx], prof_y[idx],
                                              prof_yerr[idx], unit_point)

        return cls(prof_x, prof_y, prof_yerr, **kwcopy)
