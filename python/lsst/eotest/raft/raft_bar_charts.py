"""
Code to produce bar charts for EO measurements for 9 CCDs in a raft.
"""
from __future__ import print_function, absolute_import
import os
from collections import OrderedDict
import matplotlib.pyplot as plt
import lsst.eotest.sensor as sensorTest

__all__ = ['RaftBarCharts']

class RaftBarCharts(object):
    """
    Class to produce bar charts based on measurements in eotest
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
            print('processing', os.path.basename(filename))
            self.sensor_ids[slot] = filename.split('_')[0]
            self.results[slot] = sensorTest.EOTestResults(filename)

    def make_bar_chart(self, column, ylabel=None, spec=None, step=20,
                       marker='r--', title=None):
        """
        Make a bar chart for the specified column for all 9 sensors
        in a raft.

        Parameters
        ----------
        column : str
            The name of the column to plot from the eotest results files,
            e.g., 'READ_NOISE'.
        ylabel : str, optional
            The y-axis label.  If None (default), then use the column name.
        spec : float, optional
            The value of the acceptance specification, which will be
            plotted as a line using the marker option.  If None, then
            no spec line will be plotted.
        step : float, optional
            The x-axis offset between amplifier results.  Default: 20
        marker : str, optional
            The color and line style of the plotted specification line.
            Default: 'r--'
        title : str, None
            Title of the plot. If None, then leave blank.

        Returns
        ------
        matplotlib.figure.Figure
            The figure containing the plot.
        """
        fig = plt.figure()
        for slot, results in self.results.items():
            xoffset = self._raft_slots[slot]*step
            plt.bar(results['AMP'] + xoffset, results[column], width=1)
        xtick_values = [step*i + step/2 for i in range(len(self._raft_slots))]
        plt.xticks(xtick_values, self._raft_slots.keys())
        if ylabel is None:
            ylabel = column
        plt.ylabel(ylabel)
        if spec is not None:
            xmin, xmax = plt.axis()[:2]
            plt.plot([xmin, xmax], [spec, spec], marker)
        if title is not None:
            plt.title(title)
        return fig
