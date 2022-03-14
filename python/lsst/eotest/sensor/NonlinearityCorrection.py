"""
Code to apply non-linearity correction.
"""
import copy
import warnings
import numpy as np

class NonlinearityCorrection:
    """
    Class to apply a non-linearity correction

    The point of this calls is to serve as a callable object that will
    linearize bias-subtracted data

    corrected_adu = nlc(amp, uncorrected_adu)

    This is implemented as a power-law fit to the ratio model
    signal/measured signal, where the model signal is derived from a
    linear fit (with zero y-intercept) to the measured signal vs
    incident flux data.
    """
    def __init__(self, params, fitted_range):
        """
        Parameters
        ----------
        params : dict
           Dictionary of power-law parameters to pass to np.poly1d, keyed
           by amp number.
        """
        self.params = params
        self.fitted_range = fitted_range

    def __call__(self, amp, adu):
        """Apply the non-linearity correction to the specified amp."""
        func = np.poly1d(self.params[amp])
        adu_bounds = self.fitted_range[amp]
        if not isinstance(adu, np.ndarray):
            # Consider single adu values, applying adu_bounds to
            # avoid extrapolaton.
            if adu < adu_bounds[0]:
                ratio = func(np.log(adu_bounds[0]))
            elif adu > adu_bounds[1]:
                ratio = func(np.log(adu_bounds[1]))
            else:
                ratio = func(np.log(adu))
            return ratio*adu
        # Consider numpy arrays
        adu_tmp = copy.copy(adu)
        adu_tmp[np.where(adu_tmp < adu_bounds[0])] = adu_bounds[0]
        adu_tmp[np.where(adu_tmp > adu_bounds[1])] = adu_bounds[1]
        return func(np.log(adu_tmp))*adu

    @classmethod
    def create_from_det_response(cls, detresp, gains, order=20,
                                 fit_range=(100, 3e4), adu_range=(100, 1e5),
                                 max_ratio_dev=0.02):
        """
        Create a NonlinearityCorrection object from a DetectorResponse file

        Note that the DetectorResponse files store the signal in
        electrons, but we want a correction that works on ADU, so we
        must convert using the gains.

        Parameters
        ----------
        detresp : DetectorResponse
            An object with the detector response calculated from flat-pair data.
        gains : dict
            Dictionary with gains keyed by amp number.
        order : int [20]
            Order of the polynomial to fit.
        fit_range : (float, float) [(100, 3e4)]
            Range in e- over which to fit the linear model to the
            measured signal vs flux data.
        adu_range : (float, float) [(100, 1e5)]
            Range in ADU over which to fit the polynomial model of
            the nonlinear correction.
        max_ratio_dev : float [0.02]
            Maximum deviation from unity for fitting the
            uncorrected/corrected signal level vs log(corrected signal)

        Returns
        -------
        NonlinearityCorrection
        """
        params = dict()
        fitted_range = dict()
        xdata = detresp.flux
        for amp, measured_signal in detresp.Ne.items():
            ydata = copy.copy(measured_signal)
            ypeak_index = np.argmax(ydata)
            index = np.where((fit_range[0] < ydata) & (ydata < fit_range[1])
                             & (xdata < xdata[ypeak_index]))
            # Fit linear model with zero y-intercept.
            slope = np.sqrt(sum(ydata[index]**2/xdata[index])
                            /sum(xdata[index]))
            # Fit a power-law to the ratio of model to measured
            # signals as a function of the log of the measured signal,
            # restricting to data below the location of the peak
            # signal.
            ratio = slope*xdata/ydata
            # Convert to ADU.
            ydata /= gains[amp]
            # Restrict the range of measured signal to fit and avoid
            # deviations from unity > max_ratio_dev.
            index = np.where((adu_range[0] < ydata) & (ydata < adu_range[1])
                             & (1 - max_ratio_dev < ratio)
                             & (ratio < 1 + max_ratio_dev))
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                params[amp] = np.polyfit(np.log(ydata[index]), ratio[index],
                                         order)
            fitted_range[amp] = min(ydata[index]), max(ydata[index])
        return cls(params, fitted_range)
