"""
Code to apply non-linearity correction.
"""
import copy
import numpy as np
from astropy.io import fits
import scipy.optimize
from scipy.interpolate import UnivariateSpline
from lsst.eotest.fitsTools import fitsTableFactory, fitsWriteto


def lin_func(pars, xvals):
    """Return a line whose slope is pars[0]"""
    return pars[0]*xvals

def chi2_model(pars, xvals, yvals):
    """Return the chi2 w.r.t. the model"""
    return (yvals - lin_func(pars, xvals))/np.sqrt(lin_func(pars, xvals))

def make_profile_hist(xbin_edges, xdata, ydata, yerrs=None, stderr=False):
    """Build a profile historgram

    Parameters
    ----------
    xbin_edges : `array`
        The bin edges
    xdata : `array`
        The x-axis data
    ydata : `array`
        The y-axis data
    yerrs :  `array` [None]
        The errors on the y-axis points
    stderr : `bool` [False]
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

    This is implemented as a spline interpolation for each of the
    amplifiers on a CCD
    """
    def __init__(self, prof_x, prof_y, prof_yerr, s=1e-6, ext=3):
        """
        Parameters
        ----------
        prof_x : dict
            Dictionary of nbins values for the x-axis of the correction function
        prof_y : `array`
            Dictionary of nbins values for the y-axis of the correction function
        prof_yerr : `array`
            Dictionary of nbins values for the y-axis errors of the correction
            function
        s : float [1e-6]
            Smoothing parameter for UnivariateSpline
        ext : int [3]
            Extrapolation mode of UnivariateSpline, ext=3 returns the boundary
            value.
        """
        self._prof_x = prof_x
        self._prof_y = prof_y
        self._prof_yerr = prof_yerr
        self._nxbins = len(list(prof_x.values())[0])

        self._spline_dict = {}
        for amp in prof_x:
            idx_sort = np.argsort(self._prof_x[amp])
            profile_x = self._prof_x[amp][idx_sort]
            profile_y = self._prof_y[amp][idx_sort]
            if self._prof_yerr is not None:
                profile_yerr = self._prof_yerr[amp][idx_sort]
                mask = profile_yerr >= 0.
            else:
                mask = np.ones(profile_x.shape)
            self._spline_dict[amp] = UnivariateSpline(profile_x[mask],
                                                      profile_y[mask],
                                                      s=s, ext=ext)

    def __getitem__(self, amp):
        """Get the function that corrects a particular amp"""
        return self._spline_dict[amp]

    def __call__(self, amp, adu):
        """Apply the non-linearity correction to a particular amp"""
        return adu*self._spline_dict[amp](adu)

    def write_to_fits(self, fits_file):
        """Write this object to a FITS file"""
        output = fits.HDUList()
        output.append(fits.PrimaryHDU())

        col_amps = fits.Column(name='amps', format='J', unit='None',
                               array= np.array(list(self._prof_x.keys())))

        def convert_to_array(profile_dict):
            return np.array([_ for _ in profile_dict.values()])

        col_prof_x = fits.Column(name='prof_x', format='%iE' % self._nxbins,
                                 unit='ADU', array=convert_to_array(self._prof_x))
        col_prof_y = fits.Column(name='prof_y_corr', format='%iE' % self._nxbins,
                                 unit='ADU', array=convert_to_array(self._prof_y))
        col_prof_yerr = fits.Column(name='prof_yerr', format='%iE' % self._nxbins,
                                    unit='ADU',
                                    array=convert_to_array(self._prof_yerr))

        fits_cols = [col_amps, col_prof_x, col_prof_y, col_prof_yerr]
        hdu = fitsTableFactory(fits_cols)
        hdu.name = 'nonlin'
        output.append(hdu)

        fitsWriteto(output, fits_file, overwrite=True)
        return output

    def save_plots(self, plotfile=None, ymin=None, ymax=None, figsize=(15, 10)):
        """Save plots showing the nonlinearity correction"""
        import matplotlib.pyplot as plt

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

        amp = 1
        for i_row in range(4):
            for i_col in range(4):
                axes = axs[i_row, i_col]
                if ymin is not None or ymax is not None:
                    axes.set_ylim(ymin, ymax)
                mask = self._prof_yerr[amp] >= 0.
                x_masked = self._prof_x[amp][mask]
                xline = np.linspace(1., x_masked.max(), 1001)
                model = self._spline_dict[amp](xline)
                axes.errorbar(x_masked, self._prof_y[amp][mask],
                              yerr=self._prof_yerr[amp][mask], fmt='.')
                axes.plot(xline, model, 'r-')
                amp += 1
        if plotfile is None:
            fig.show()
        else:
            fig.savefig(plotfile)

    @classmethod
    def create_from_table(cls, table, s=1e-6, ext=3):
        """Create a NonlinearityCorrection object from a fits file

        Parameters
        ----------
        table : `Table`
            The table data used to build the nonlinearity correction
        s : float [1e-6]
            Smoothing parameter for UnivariateSpline
        ext : int [3]
            Extrapolation mode of UniveriateSpline. ext=3 returns the boundary
            value

        Returns
        -------
        nl : `NonlinearityCorrection`
            The requested object
        """
        amps = table.data['amps']
        prof_x = dict(zip(amps, table.data['prof_x']))
        prof_y = dict(zip(amps, table.data['prof_y_corr']))
        prof_yerr = dict(zip(amps, table.data['prof_yerr']))
        return cls(prof_x, prof_y, prof_yerr, s=s, ext=ext)

    @classmethod
    def create_from_fits_file(cls, fits_file, hdu_name='nonlin', s=1e-6,
                              ext=3):
        """Create a NonlinearityCorrection object from a fits file

        Parameters
        ----------
        fits_file : `str`
            The file with the data used to build the nonlinearity correction
        hdu_name : `str` ['nonlin']
            The name of the HDU with the nonlinearity correction data
        s : float [1e-6]
            Smoothing parameter for UnivariateSpline
        ext : int [3]
            Extrapolation mode of UniveriateSpline. ext=3 returns the boundary
            value

        Returns
        -------
        nl : `NonlinearityCorrection`
            The requested object
        """
        hdulist = fits.open(fits_file)
        table = hdulist[hdu_name]
        nl = cls.create_from_table(table, s=s, ext=ext)
        hdulist.close()
        return nl

    @staticmethod
    def _correct_fixed_point(profile_x, profile_y, profile_yerr, fixed_point):
        """Force the spline to go through one at a particular x-value

        Parameters
        ----------
        profile_x : `array`
            The x-bin centers
        profile_y : `array`
            The b-bin values
        profile_yerr : `array`
            The y-bin errors
        fixed_point : `float`
            The x-value where the spline should go through one

        Returns
        -------
        y_vals_corr
            The adjusted y-values
        y_errs_corr
            The adjusted y-errors
        """
        uni_spline = UnivariateSpline(profile_x, profile_y)
        y_value_0 = uni_spline(fixed_point)

        y_vals_corr = profile_y/y_value_0
        y_errs_corr = profile_yerr
        return y_vals_corr, y_errs_corr

    @classmethod
    def create_from_det_response(cls, detresp, gains, fit_range=(0, 9e4),
                                 nprofile_bins=10, fixed_point=None,
                                 max_ratio_dev=0.02, s=1e-6, ext=3):
        """Create a NonlinearityCorrection object DetectorResponse FITS file

        Note that the DetectorResponse files typically store the
        signal in electrons, but we want a correction that works on
        ADU, so we have to remove the gains.

        Parameters
        ----------
        detresp : `DetectorResponse`
            An object with the detector response calculated from flat-pair files
        gains : `dict`
            Dictionary with amp gains keyed by amp number.
        fit_range : `tuple` [(0., 9e4)]
            The range over which to define the non-linearity
        nprofile_bins : `int` [10]
             The number of bins to use in the profile.
        fixed_point : `float` [None]
             X-value at which the correction should vanish, defaults to 0.
             If `None` then this will simply use the pivot point of the fit
             to the data
        max_ratio_dev : float [0.02]
             Maximum deviation from unity of ratio points to include in
             spline fit.
        s : float [1e-6]
            Smoothing parameter for UnivariateSpline
        ext : int [3]
            Extrapolation mode of UniveriateSpline. ext=3 returns the boundary
            value

        Returns
        -------
        nl : `NonlinearityCorrection`
            The requested object

        """
        prof_x = dict()
        prof_y = dict()
        prof_yerr = dict()

        # For each amp, fit a linear model to the signal vs
        # incident flux.
        xdata = detresp.flux
        for amp in detresp.Ne:
            ydata = copy.copy(detresp.Ne[amp])
            # The nominal fit_range applies to e-/pixel.  We also want
            # to avoid fitting data past the saturation peak.
            ypeak_index = np.argmax(ydata)
            x_at_ypeak = xdata[ypeak_index]
            mask = np.where((fit_range[0] < ydata) & (fit_range[1] > ydata)
                            & (xdata < x_at_ypeak))
            # Convert back to ADU
            ydata /= gains[amp]

            xdata_fit = xdata[mask]
            ydata_fit = ydata[mask]
            mean_slope = (ydata_fit/xdata_fit).mean()
            pars = (mean_slope,)
            results = scipy.optimize.leastsq(chi2_model, pars,
                                             full_output=1,
                                             args=(xdata_fit, ydata_fit))

            # Compute the ratio of the linear model to the measured
            # signal in ADU.  This ratio would be the correction
            # factor assuming the y-value computed from the linear fit
            # is the desired signal.  The correction factor would then
            # be a function of the measured signal in ydata.
            ratio = lin_func(results[0], xdata)/ydata

            # Avoid ratio points > +/-max_ratio_dev from unity and
            # don't try to fit a spline past the saturation peak.
            index = np.where((1 - max_ratio_dev < ratio)
                             & (ratio < 1 + max_ratio_dev)
                             & (xdata < x_at_ypeak))
            ydata = ydata[index]
            ratio = ratio[index]

            # Bin the data over the full range of signal values that
            # pass the above cuts.
            xbins = np.linspace(min(ydata), max(ydata), nprofile_bins + 1)
            prof_x[amp], prof_y[amp], prof_yerr[amp] \
                = make_profile_hist(xbins, ydata, ratio, stderr=True)

            if fixed_point is not None:
                prof_y[amp], prof_yerr[amp] \
                    = cls._correct_fixed_point(prof_x[amp], prof_y[amp],
                                               prof_yerr[amp], fixed_point)

        return cls(prof_x, prof_y, prof_yerr, s=s, ext=ext)
