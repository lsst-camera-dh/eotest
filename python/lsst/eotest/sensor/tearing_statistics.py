"""
Code to compute Pierre Antilogus' tearing detection statistics.
See https://jira.slac.stanford.edu/browse/LSSTTD-1273.
"""

from __future__ import print_function
import os
from collections import namedtuple, OrderedDict
import numpy as np
import matplotlib.pyplot as plt
from .MaskedCCD import MaskedCCD

__all__ = ['TearingStats', 'AmpTearingStats']

class TearingStats(OrderedDict):
    """
    Tearing statistics for a CCD, keyed by amp number, 1-17.
    """
    def __init__(self, fitsfile, buf=10):
        """
        Parameters
        ----------
        fitsfile: str
            Single sensor FITS file.
        buf: int [10]
            Number of pixels to avoid on leading and trailing edge of
            serial overscan to compute the bias level for each row.
        """
        super(TearingStats, self).__init__()
        self.fitsfile = fitsfile
        ccd = MaskedCCD(fitsfile)
        for amp in ccd:
            self[amp] = AmpTearingStats(ccd[amp], ccd.amp_geom, buf=buf)

    def has_tearing(self, cut1=0.05, cut2=-0.01, nsig=1, ntear_min=10):
        """
        Returns True if tearing is found given the ratio cuts, signal
        level, and mininum number of individual detections over all
        amplifier edges.

        Parameters
        ----------
        cut1: float [0.05]
        cut2: float [-0.01]
        nsig: float [1]
        ntear_min: int [10]
        """
        ntear = 0
        for amp in self:
            rstats1, rstats2 = self[amp].rstats
            if (rstats1.diff - cut1 > nsig*rstats1.error and
                rstats2.diff - cut2 > nsig*rstats2.error):
                ntear += 1
            if (rstats1.diff - cut2 > nsig*rstats1.error and
                rstats2.diff - cut1 > nsig*rstats2.error):
                ntear += 1
        return ntear > ntear_min

    def plot_profiles(self, fig=None, figsize=(10, 10), title=None):
        """
        Plot the tearing profiles for both edges of each amplifier.

        Parameters
        ----------
        fig: matplotlib.Figure [None]
             If None, then one will be created.
        figsize: tuple(int, int) [10, 10]
             Size of the figure in inches.
        title: str [None]
             Title of the figure. If None, then the basename of the FITS
             file will be used along with the text of "tearing" or
             "no tearing" based on the analysis outcome.
        """
        saved_figsize = plt.rcParams['figure.figsize']
        saved_fontsize = plt.rcParams['font.size']
        try:
            if fig is None:
                plt.rcParams['figure.figsize'] = figsize
                plt.rcParams['font.size'] = 10
                fig = plt.figure()
            frame_axes = fig.add_subplot(111, frameon=False)
            if title is None:
                title = os.path.basename(self.fitsfile)
            outcome = "tearing detected" if self.has_tearing() else "no tearing"
            frame_axes.set_title(title + ": " + outcome)
            frame_axes.set_xlabel('\ny-pixel', fontsize=12)
            frame_axes.set_ylabel('ratio of counts for outer two columns\n\n',
                                  fontsize=12)
            frame_axes.get_xaxis().set_ticks([])
            frame_axes.get_yaxis().set_ticks([])
            for amp in self:
                prof1, prof2 = self[amp].ratio_profiles
                ratios1, ratios2 = self[amp].ratios
                stdevs1, stdevs2 = self[amp].stdevs
                ax = fig.add_subplot(4, 4, amp)
                plt.errorbar(range(len(prof1)), prof1, fmt='.', color='green',
                             alpha=0.3, label='first edge', zorder=1)
                plt.errorbar(range(len(prof2)), prof2, fmt='.', color='blue',
                             alpha=0.3, label='last edge', zorder=1)
                plt.errorbar(self[amp].ylocs, ratios1, yerr=stdevs1,
                             xerr=3*[self[amp].dy], fmt='.', color='black',
                             zorder=10, markersize=1)
                plt.errorbar(self[amp].ylocs, ratios2, yerr=stdevs2,
                             xerr=3*[self[amp].dy], fmt='.', color='black',
                             zorder=10, markersize=1)
                ax.annotate('amp %d' % amp, (0.65, 0.9),
                            xycoords='axes fraction', fontsize='small')
            plt.tight_layout()
        finally:
            plt.rcParams['figure.figsize'] = saved_figsize
            plt.rcParams['font.size'] = saved_fontsize


class AmpTearingStats(object):
    """
    Class to compute Pierre Antilogus' tearing statistics based on
    the size of the transitions in the ratio of the two outermost
    columns of pixels on both sides of an amplifier.
    """
    def __init__(self, full_segment, amp_geom, buf=10):
        """
        Parameters
        ----------
        full_segment: lsst.afw.image.Image
            Image object of the full segment of an amplifier.
        amp_geom: lsst.eotest.sensor.AmplifierGeometry
            Object containing the amplifier pixel geometry.
        buf: int [10]
            Number of pixels to avoid on leading and trailing edge of
            serial overscan to compute the bias level for each row.
        """
        self._subtract_overscan(full_segment.getImage(), amp_geom, buf)
        self._ratio_profiles = None
        self._ratios = None
        self._stdevs = None
        self._rstats = None
        self.ylocs = 150, 1000, 1950
        self.dy = 50

    def _subtract_overscan(self, full_segment, amp_geom, buf):
        oscan = full_segment.Factory(full_segment, amp_geom.serial_overscan)
        imarr = full_segment.Factory(full_segment, amp_geom.imaging).getArray()
        ny, nx = imarr.shape
        bias = np.array(nx*[oscan.getArray()[:ny, buf:-buf].mean(axis=1)])
        self.imarr = imarr - bias.transpose()

    @property
    def rstats(self):
        """Full range of mean edge ratio values with errors."""
        RatioStats = namedtuple('RatioStats', ('diff', 'error'))
        if self._rstats is None:
            self._rstats = []
            for i in range(2):
                r = self.ratios[i]
                dr = self.stdevs[i]
                jmin = np.argmin(r)
                jmax = np.argmax(r)
                self._rstats.append(
                    RatioStats(r[jmax] - r[jmin],
                               np.sqrt(dr[jmax]**2 + dr[jmin]**2)/10.))
        return self._rstats

    @property
    def ratio_profiles(self):
        """Profiles of the edge pixel ratios."""
        if self._ratio_profiles is None:
            self._ratio_profiles = (self.imarr[:, 1]/self.imarr[:, 0],
                                    self.imarr[:, -2]/self.imarr[:, -1])
        return self._ratio_profiles

    @property
    def ratios(self):
        """Mean values of ratio profiles evaluated at self.yloc locations."""
        if self._ratios is None:
            self._compute_ratios()
        return self._ratios

    @property
    def stdevs(self):
        """Standard deviatoins of ratio profiles evaluated at self.yloc
        locations."""
        if self._stdevs is None:
            self._compute_ratios()
        return self._stdevs

    def _compute_ratios(self):
        self._ratios = []
        self._stdevs = []
        for profile in self.ratio_profiles:
            my_ratios = []
            my_stdevs = []
            for yloc in self.ylocs:
                ymin, ymax = yloc - self.dy, yloc + self.dy
                my_ratios.append(profile[ymin:ymax].mean())
                my_stdevs.append(profile[ymin:ymax].std())
            self._ratios.append(my_ratios)
            self._stdevs.append(my_stdevs)

if __name__ == '__main__':
    plt.ion()

    fitsfile = '00_FlatBaseline_0000_20180426205147.fits.gz'
    ts = TearingStats(fitsfile)
    ts.plot_profiles()
    print("has tearing:", ts.has_tearing())
