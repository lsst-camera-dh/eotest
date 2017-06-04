"""
Code to make a summary of QE performance at the raft-level.
"""
from collections import OrderedDict
import json
import itertools
import astropy.io.fits as fits
import matplotlib.pyplot as plt

__all__ = ['qe_summary_plot', 'parse_summary_lims', 'parse_qe_files']

def parse_qe_files(qe_files):
    """
    Extract the QE results from the results files produced by the
    single sensore QeTask.

    Parameters
    ----------
    qe_files : dict
        Dictionary, keyed by slot, of the filenames of the QE results files.

    Returns
    -------
    OrderedDict(dict) : ordered dict keyed, by band
        (ugrizy), of ordered dicts, keyed by slot, of the QE values
        for that sensor in percentages.
    """
    slots = ['S%i%i' % pair for pair in itertools.product(range(3), range(3))]

    qe = OrderedDict()
    for band in 'ugrizy':
        qe[band] = OrderedDict([(slot, 0) for slot in slots])
    for slot, qe_file in qe_files.items():
        qe_data = fits.open(qe_file)['QE_BANDS'].data
        for band, value in zip(qe_data['BAND'], qe_data['DEVICE_MEAN']):
            qe[band][slot] = value
    return qe

def parse_summary_lims(summary_lims):
    """
    Parse the summary.lims file that contains the raft-level QE results.

    Parameters
    ----------
    summary_lims : str
        summary.lims file containing the raft-level QE results for
        each sensor.

    Returns
    -------
    OrderedDict(dict) : ordered dict keyed, by band
        (ugrizy), of ordered dicts, keyed by slot, of the QE values
        for that sensor in percentages.
    """
    slots = ['S%i%i' % pair for pair in itertools.product(range(3), range(3))]

    qe = OrderedDict()
    for band in 'ugrizy':
        qe[band] = OrderedDict([(slot, 0) for slot in slots])

    qe_results = json.load(open(summary_lims))
    for item in qe_results:
        if item.has_key('QE'):
            qe[item['band']][item['slot']] = item['QE']

    return qe

def qe_summary_plot(qe, figsize=(6, 8), title='', qe_spec=None):
    """Function to make a summary plot of QE values in each band for each
    of sensors in a raft.

    Parameters
    ----------
    qe : OrderedDict
        Ordered dict keyed, by band (ugrizy), of dictionaries, keyed
        by slot, of the QE values for that sensor in percentages.
    figsize : tuple, optional
        Figure size in inches.  Default: (6, 8)
    title : str, optional
        Overall figure title.  Default: ''
    qe_spec : dict, optional
        Dictionary of throughput specification values, keyed by 'ugrizy'.
        By default the LCA-57 values are used.

    """
    if qe_spec is None:
        # Set spec values to LCA-57 values.
        qe_spec = dict(u=41, g=78, r=83, i=82, z=75, y=20.1)

    plt.rcParams['figure.figsize'] = figsize
    fig = plt.figure()
    frame_axes = fig.add_subplot(111, frameon=False)
    frame_axes.set_title(title)

    frame_axes.set_xlabel('\nslot')
    frame_axes.get_xaxis().set_ticks([])
    frame_axes.get_yaxis().set_ticks([])

    bounds = -1, 9, 0, 100
    slots = qe.values()[0].keys()
    xtick_values = range(len(slots))
    for i, qe_vals in enumerate(qe.items()):
        band = qe_vals[0]
        fig.add_subplot(6, 1, i+1)
        plt.errorbar(xtick_values, qe_vals[1].values(), fmt='.')
        plt.xticks(xtick_values, len(slots)*[''])
        plt.axis(bounds)
        plt.ylabel('%s QE' % band)
        plt.plot(bounds[:2], [qe_spec[band], qe_spec[band]], 'r--')
        plt.grid()
    plt.xticks(xtick_values, slots)
    return fig
