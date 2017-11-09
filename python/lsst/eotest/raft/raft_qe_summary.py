"""
Code to make a summary of QE performance at the raft-level.
"""
from builtins import range
from collections import OrderedDict
import json
import itertools
import matplotlib.pyplot as plt

__all__ = ['qe_summary_plot']


def qe_summary_plot(summary_lims, figsize=(6, 8), title='', qe_spec=None):
    """
    Function to make a summary plot of QE values in each band for each
    of sensors in a raft.

    Parameters
    ----------
    summary_lims : str
        summary.lims file containing the raft-level QE results for
        each sensor.
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

    qe_results = json.load(open(summary_lims))

    slots = ['S%i%i' % pair for pair in itertools.product(list(range(3)), list(range(3)))]

    qe = OrderedDict()
    for band in 'ugrizy':
        qe[band] = OrderedDict([(slot, 0) for slot in slots])

    for item in qe_results:
        if 'QE' in item:
            qe[item['band']][item['slot']] = item['QE']

    plt.rcParams['figure.figsize'] = figsize
    fig = plt.figure()
    frame_axes = fig.add_subplot(111, frameon=False)
    frame_axes.set_title(title)

    frame_axes.set_xlabel('\nslot')
    frame_axes.get_xaxis().set_ticks([])
    frame_axes.get_yaxis().set_ticks([])

    bounds = -1, 9, 0, 100
    xtick_values = list(range(len(slots)))
    for i, qe_vals in enumerate(qe.items()):
        band = qe_vals[0]
        fig.add_subplot(6, 1, i+1)
        plt.errorbar(xtick_values, list(qe_vals[1].values()), fmt='.')
        plt.xticks(xtick_values, len(slots)*[''])
        plt.axis(bounds)
        plt.ylabel('%s QE' % band)
        plt.plot(bounds[:2], [qe_spec[band], qe_spec[band]], 'r--')
        plt.grid()
    plt.xticks(xtick_values, slots)
    return fig
