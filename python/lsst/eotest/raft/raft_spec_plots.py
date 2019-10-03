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
        for slot, filename in list(results_files.items()):
            self.sensor_ids[slot] = filename.split('_')[0]
            self.results[slot] = sensorTest.EOTestResults(filename)

    def _draw_slot_boundary(self, slot, step=20, namps=16, marker='k:'):
        if slot == 'S00':
            plt.axvline(-2, linestyle=marker[1:], color=marker[0])
        xbound = (step - namps)/2. + step*self._raft_slots[slot] + namps + 1
        plt.axvline(xbound, linestyle=marker[1:], color=marker[0])

    @staticmethod
    def _apply_ybounds(yvalues, ybounds):
        if ybounds is None:
            return yvalues
        my_yvalues = np.zeros(len(yvalues))
        for i, yy in enumerate(yvalues):
            my_yvalues[i] = min(max(yy, ybounds[0]), ybounds[1])
        return my_yvalues

    def make_plot(self, column, ylabel=None, spec=None, step=20, yscaling=1,
                  marker='r--', title=None, ylog=False, figsize=(8, 6),
                  ybounds=None):
        """
        Make a plot for the specified column for all 9 sensors
        in a raft.

        Parameters
        ----------
        column : str
            The name of the column to plot from the eotest results files,
            e.g., 'READ_NOISE'.
        ylabel : str [None]
            The y-axis label.  If None, then use the column name.
        spec : float or sequence [None]
            The value(s) of the acceptance specification(s), which will be
            plotted as a line using the marker option.  If None, then
            no spec line will be plotted.
        step : float [20]
            The x-axis offset between amplifier results.
        yscaling : float [1]
            Scale factor to apply to y-values.
        marker : str ['r--']
            The color and line style of the plotted specification line.
        title : str [None]
            Title of the plot. If None, then leave blank.
        ylog : bool [False]
            If True, make the y-axis log scale.
        figsize : tuple(float, float) [(8, 6)]
            Figure size in inches.
        ybounds : tuple(float, float) [None]
            Outer range of y-axis bounds.  If None, then default matplotlib
            autoscaling is used; otherwise, these bounds are only enforced
            if the corresponding autoscale bound exceeds this outer range.

        Returns
        ------
        matplotlib.figure.Figure
            The figure containing the plot.
        """
        plt.rcParams['figure.figsize'] = figsize
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        if ybounds is not None:
            ybounds = yscaling*ybounds[0], yscaling*ybounds[1]
        for slot, results in list(self.results.items()):
            xoffset = self._raft_slots[slot]*step
            yvalues = self._apply_ybounds(yscaling*results[column], ybounds)
            plt.plot(results['AMP'] + xoffset, yvalues, 'b.')
        xtick_values = [step*i + step/2 for i in range(len(self._raft_slots))]
        plt.xticks(xtick_values, list(self._raft_slots.keys()))
        if ylabel is None:
            ylabel = column
        plt.ylabel(ylabel)
        if spec is not None:
            if not hasattr(spec, '__iter__'):
                spec = (spec,)
            for spec_value in spec:
                yval = yscaling*spec_value
                plt.axhline(yval, color=marker[0], linestyle=marker[1:])
        if title is not None:
            plt.title(title)
        if ylog:
            ax.set_yscale('log', nonposy='clip')
        if ybounds is not None:
            axis = list(plt.axis())
            axis[-2] = max(axis[-2], ybounds[0])
            axis[-1] = min(axis[-1], ybounds[1])
            plt.axis(axis)
        for slot in self.results:
            self._draw_slot_boundary(slot, step=step)
        return fig

    def make_multi_column_plot(self, columns, ylabel=None, spec=None, step=20,
                               yscaling=1, yerrors=False, marker='r--',
                               title=None, add_legend=True, ylog=False,
                               colors=None, figsize=(8, 6), ybounds=None):
        """
        Make a plot for the specified columns for all 9 sensors
        in a raft.

        Parameters
        ----------
        column : tuple(str)
            The names of the columns to plot from the eotest results files,
            e.g., ('GAIN', 'PTC_GAIN')
        ylabel : str [None]
            The y-axis label.  If None, then use the column name.
        spec : float or sequence [None]
            The value(s) of the acceptance specification(s), which will be
            plotted as a line using the marker option.  If None, then
            no spec line will be plotted.
        step : float [20]
            The x-axis offset between amplifier results.
        yscaling : float [1]
            Scale factor to apply to y-values.
        yerrors : bool [False]
            Flag to plot errors, assuming the column name is of the form
            <column>'_ERROR', e.g., 'GAIN_ERROR'.  Default: False
        marker : str ['r--']
            The color and line style of the plotted specification line.
        title : str [None]
            Title of the plot. If None, then leave blank.
        add_legend : bool [True]
            If True, then add a legend.
        ylog : bool [False]
            If True, make the y-axis log scale.
        colors : tuple [None]
            Colors to cycle over, e.g., 'krgbcym'.  If None, then use
            the default color cycle.
        figsize : tuple(float, float) [(8, 6)]
            Figure size in inches.
        ybounds : tuple(float, float) [None]
            Outer range of y-axis bounds.  If None, then default matplotlib
            autoscaling is used; otherwise, these bounds are only enforced
            if the corresponding autoscale bound exceeds this outer range.

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
        if ybounds is not None:
            ybounds = yscaling*ybounds[0], yscaling*ybounds[1]
        for icol, column in enumerate(columns):
            x, y = [], []
            yerr = []
            for slot, results in list(self.results.items()):
                xoffset = self._raft_slots[slot]*step
                x.extend(results['AMP'] + xoffset + icol*dx)
                y.extend(yscaling*results[column])
                if yerrors:
                    yerr.extend(yscaling*results[column + '_ERROR'])
            color = next(color_cycler)['color']
            y = self._apply_ybounds(y, ybounds)
            plt.plot(x, y, '.', label=column, color=color)
            if yerrors:
                plt.errorbar(x, y, fmt='.', yerr=yerr, color=color)
        xtick_values = [step*i + step/2 for i in range(len(self._raft_slots))]
        plt.xticks(xtick_values, list(self._raft_slots.keys()))
        if ylabel is None:
            ylabel = column
        plt.ylabel(ylabel)
        if spec is not None:
            if not hasattr(spec, '__iter__'):
                spec = (spec,)
            for spec_value in spec:
                yval = yscaling*spec_value
                plt.axhline(yval, color=marker[0], linestyle=marker[1:])
        if title is not None:
            plt.title(title)
        if add_legend:
            plt.legend(loc=0)
        if ylog:
            ax.set_yscale('log', nonposy='clip')
        if ybounds is not None:
            axis = list(plt.axis())
            axis[-2] = max(axis[-2], ybounds[0])
            axis[-1] = min(axis[-1], ybounds[1])
            plt.axis(axis)
        for slot in self.results:
            self._draw_slot_boundary(slot, step=step)
        self._draw_slot_boundary('S00', step=step)
        return fig
