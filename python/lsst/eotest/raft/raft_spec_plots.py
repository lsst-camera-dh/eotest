"""
Code to produce plots of EO measurements for 9 CCDs in a raft.
"""
from __future__ import print_function, absolute_import
import os
from collections import OrderedDict
import cycler
import numpy as np
import matplotlib.pyplot as plt
import lsst.eotest.sensor as sensorTest

__all__ = ['RaftSpecPlots']

class RaftSpecPlots(object):
    """
    Class to produce plots of measured specifications in eotest
    results files.
    """
    _raft_slots = \
        OrderedDict([(slot, i) for i, slot in
                     enumerate('S00 S01 S02 S10 S11 S12 S20 S21 S22'.split())])
    def __init__(self, results_files):
        """
        Constructor.

        Parameters
        ----------
        results_files : dict
        """
        self.results = dict()
        self.sensor_ids = dict()
        for slot, filename in results_files.items():
            self.sensor_ids[slot] = filename.split('_')[0]
            self.results[slot] = sensorTest.EOTestResults(filename)

    def _draw_slot_boundary(self, slot, step=20, namps=16, marker='k:'):
        xbound = (step - namps)/2. + step*self._raft_slots[slot] + namps + 1
        axis = plt.axis()
        ymin, ymax = axis[2:]
        plt.plot([xbound, xbound], [ymin, ymax], marker)
        plt.axis(axis)

    def make_plot(self, column, ylabel=None, spec=None, step=20, yscaling=1,
                  marker='r--', title=None, ylog=False, figsize=(8, 6),
                  ymax=None, yerrors=False):
        """
        Make a plot for the specified column for all 9 sensors
        in a raft.

        Parameters
        ----------
        column : str
            The name of the column to plot from the eotest results files,
            e.g., 'READ_NOISE'.
        ylabel : str, optional
            The y-axis label.  If None (default), then use the column name.
        spec : float or sequence, optional
            The value(s) of the acceptance specification(s), which will be
            plotted as a line using the marker option.  If None, then
            no spec line will be plotted.
        step : float, optional
            The x-axis offset between amplifier results.  Default: 20
        yscaling : float, optional
            Scale factor to apply to y-values.  Default=1.
        marker : str, optional
            The color and line style of the plotted specification line.
            Default: 'r--'
        title : str, optional
            Title of the plot. If None, then leave blank.
        ylog : bool, optional
            If True, make the y-axis log scale. Default: False
        figsize : tuple(float, float), optional
            Figure size in inches. Default: (8, 6)
        ymax : float, optional
            Upper bound of plot in data coordinates.  If not None (default),
            then use the max(plt.axis()[-1], ymax) as the upper bound
            of the plotting area.
        yerrors : bool, optional
            Flag to plot errors, assuming the column name is of the form
            <column>'_ERROR', e.g., 'GAIN_ERROR'.  Default: False

        Returns
        ------
        matplotlib.figure.Figure
            The figure containing the plot.
        """
        plt.rcParams['figure.figsize'] = figsize
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        for slot, results in self.results.items():
            xoffset = self._raft_slots[slot]*step
            x = results['AMP'] + xoffset
            y = yscaling*results[column]
            plt.plot(x, y, 'b.')
            if yerrors:
                yerr = yscaling*results[column+'_ERROR']
                plt.errorbar(x, y, fmt='.', yerr=yerr, color='blue')
        xtick_values = [step*i + step/2 for i in range(len(self._raft_slots))]
        plt.xticks(xtick_values, self._raft_slots.keys())
        if ylabel is None:
            ylabel = column
        plt.ylabel(ylabel)
        if spec is not None:
            if not hasattr(spec, '__iter__'):
                spec = (spec,)
            for spec_value in spec:
                xmin, xmax = plt.axis()[:2]
                yval = yscaling*spec_value
                plt.plot([xmin, xmax], [yval, yval], marker)
        if title is not None:
            plt.title(title)
        if ylog:
            ax.set_yscale('log', nonposy='clip')
        if ymax is not None:
            axis = list(plt.axis())
            axis[-1] = max(axis[-1], yscaling*ymax)
            plt.axis(axis)
        for slot in self.results:
            self._draw_slot_boundary(slot, step=step)
        return fig

    def make_multi_column_plot(self, columns, ylabel=None, spec=None, step=20,
                               yscaling=1, yerrors=False, marker='r--',
                               title=None, add_legend=True, ylog=False,
                               colors=None, figsize=(8, 6), ymax=None):
        """
        Make a plot for the specified columns for all 9 sensors
        in a raft.

        Parameters
        ----------
        column : tuple(str)
            The names of the columns to plot from the eotest results files,
            e.g., ('GAIN', 'PTC_GAIN')
        ylabel : str, optional
            The y-axis label.  If None (default), then use the column name.
        spec : float or sequence, optional
            The value(s) of the acceptance specification(s), which will be
            plotted as a line using the marker option.  If None, then
            no spec line will be plotted.
        step : float, optional
            The x-axis offset between amplifier results.  Default: 20
        yscaling : float, optional
            Scale factor to apply to y-values.  Default: 1.
        yerrors : bool, optional
            Flag to plot errors, assuming the column name is of the form
            <column>'_ERROR', e.g., 'GAIN_ERROR'.  Default: False
        marker : str, optional
            The color and line style of the plotted specification line.
            Default: 'r--'
        title : str, optional
            Title of the plot. If None, then leave blank.
        add_legend : bool, optional
            If True (default), then add a legend.
        ylog : bool, optional
            If True, make the y-axis log scale. Default: False
        colors : tuple, optional
            Colors to cycle over, e.g., 'krgbcym'.  If None (default),
            then use the default color cycle.
        figsize : tuple(float, float), optional
            Figure size in inches. Default: (8, 6)
        ymax : float, optional
            Upper bound of plot in data coordinates.  If not None (default),
            then use the max(plt.axis()[-1], ymax) as the upper bound
            of the plotting area.

        Returns
        ------
        matplotlib.figure.Figure
            The figure containing the plot.
        """
        plt.rcParams['figure.figsize'] = figsize
        dx = 1./len(columns)
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        # Set up the color cycling.
        if colors is not None:
            color_cycler = cycler.cycler('color', colors)()
        else:
            color_cycler = plt.rcParams['axes.prop_cycle']()
        for icol, column in enumerate(columns):
            x, y = [], []
            yerr = []
            for slot, results in self.results.items():
                xoffset = self._raft_slots[slot]*step
                x.extend(results['AMP'] + xoffset + icol*dx)
                y.extend(yscaling*results[column])
                if yerrors:
                    yerr.extend(yscaling*results[column + '_ERROR'])
            color = next(color_cycler)['color']
            plt.plot(x, y, '.', label=column, color=color)
            if yerrors:
                plt.errorbar(x, y, fmt='.', yerr=yerr, color=color)
        xtick_values = [step*i + step/2 for i in range(len(self._raft_slots))]
        plt.xticks(xtick_values, self._raft_slots.keys())
        if ylabel is None:
            ylabel = column
        plt.ylabel(ylabel)
        if spec is not None:
            if not hasattr(spec, '__iter__'):
                spec = (spec,)
            for spec_value in spec:
                xmin, xmax = plt.axis()[:2]
                yval = yscaling*spec_value
                plt.plot([xmin, xmax], [yval, yval], marker)
        if title is not None:
            plt.title(title)
        if add_legend:
            plt.legend(loc=0)
        if ylog:
            ax.set_yscale('log', nonposy='clip')
        if ymax is not None:
            axis = list(plt.axis())
            axis[-1] = max(axis[-1], yscaling*ymax)
            plt.axis(axis)
        for slot in self.results:
            self._draw_slot_boundary(slot, step=step)
        self._draw_slot_boundary('S00', step=step)
        return fig
