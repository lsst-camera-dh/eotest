"""
The gain_sequence function fits the gain for each amp as a
function of frame sequence number to the DN distributions produced by
the Fe55Task.
"""
import warnings
import numpy as np
import pandas as pd
from astropy.io import fits
from .Fe55GainFitter import Fe55GainFitter

__all__ = ['gain_sequence']


def gain_sequence(det_name, psf_results_file, chiprob_min=0.1, logger=None):
    """
    Compute the gain for each amp from the fitted DN values as a
    function of sequence number.

    Parameters
    ----------
    det_name: str
        Detector name, e.g., 'R22_S11', or sensor LSST serial number.
        This will be a column in the output data frame so that data
        frames from different sensors can be appended.
    psf_results_file: str
        PSF results file generated for each sensor by the Fe55Task.
    chiprob_min: float [0.1]
        Mininum chi-square probability to apply to fitted DN values.
    logger: logging.Logger [None]
        Logging object which will print the sequence number as data
        from each frame is fit.

    Returns
    -------
    pandas.DataFrame containing columns det_name, amp, gain, seqnum.
    """
    warnings.simplefilter('ignore')
    with fits.open(psf_results_file) as cluster_data:
        seqnums = sorted(list(set(cluster_data[1].data['SEQNUM'])))
        all_amps = range(1, cluster_data[0].header['NAMPS'] + 1)
        det_names, amps, gains, seq_nums = [], [], [], []
        for seqnum in seqnums:
            if logger is not None:
                logger.info(f'{seqnum}')
            for amp in all_amps:
                chiprob = cluster_data[amp].data['CHIPROB']
                index = np.where((chiprob > chiprob_min) &
                                 (cluster_data[amp].data['SEQNUM'] == seqnum))
                dn = cluster_data[amp].data['DN'][index]
                fitter = Fe55GainFitter(dn)
                fitter.fit()
                det_names.append(det_name)
                amps.append(amp)
                gains.append(fitter.gain)
                seq_nums.append(seqnum)
    return pd.DataFrame({'det_name': det_names, 'amp': amps,
                         'gain': gains, 'seqnum': seq_nums})
