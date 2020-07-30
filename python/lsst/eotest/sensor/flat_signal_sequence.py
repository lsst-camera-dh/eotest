"""
The flat_signal_sequence function computes the median counts per amp as
function of frame sequence number.
"""
import warnings
from collections import defaultdict
import numpy as np
import pandas as pd
import lsst.afw.math as afw_math
from .MaskedCCD import MaskedCCD
from .flatPairTask import mondiode_value

__all__ = ['flat_signal_sequence']

def flat_signal_sequence(flat_files, bias_frame=None, dark_frame=None,
                         mask_files=(), mondiode_func=mondiode_value,
                         verbose=True):
    """
    The flat_signal_sequence function computes the median counts per
    amp as function of frame sequence number.

    Parameters
    ----------
    flat_files: list-like
        List of flat FITS files to process.
    bias_frame: str [None]
        Filename of bias frame to use for bias subtraction.
    dark_frame: str [None]
        Filename of dark frame to use for dark subtraction.
    mask_files: list-like [()]
        List of mask files to use.
    mondiode_func: function object [lsst.eotest.sensor.mondiode_value]
        Function to use for computing the monitoring diode current.
    verbose: bool [True]
        If True, print the name of the file being processed

    Returns
    -------
    pandas.DataFrame containing arrays of mjd, test seqnum, integrated flux,
        and median ADU values for each amplifier.
    """
    flux_dict = dict()
    mjd_dict = dict()
    adus = None
    for item in flat_files:
        if verbose:
            print(item)
        ccd = MaskedCCD(item, bias_frame=bias_frame, dark_frame=dark_frame,
                        mask_files=mask_files)
        if adus is None:
            adus = defaultdict(lambda: {amp: 0 for amp in ccd})
        seqnum = ccd.md.get('TSEQNUM')
        mjd_dict[seqnum] = ccd.md.get('MJD-OBS')
        exptime = ccd.md.get('EXPTIME')
        flux_dict[seqnum] = mondiode_func(item, exptime)*exptime
        for amp in ccd:
            image = ccd.unbiased_and_trimmed_image(amp)
            stats = afw_math.makeStatistics(image, afw_math.MEDIAN,
                                            ccd.stat_ctrl)
            adus[seqnum][amp] = stats.getValue(afw_math.MEDIAN)

    data = {'mjd': [mjd_dict[_] for _ in adus.keys()],
            'seqnum': list(adus.keys()),
            'flux': [flux_dict[_] for _ in adus.keys()]}

    for amp in ccd:
        data[f'amp{amp:02d}'] = [_[amp] for _ in adus.values()]

    return pd.DataFrame(data=data)
