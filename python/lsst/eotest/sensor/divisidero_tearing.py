"""
Code to perform Divisadero tearing analysis.  This is slightly
revised code originally from Aaron Roodman. See LSSTTD-1440.
"""
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from astropy import stats
import lsst.eotest.image_utils as imutils
from .MaskedCCD import MaskedCCD
from .AmplifierGeometry import makeAmplifierGeometry


__all__ = ['ana_divisidero_tearing']


def normed_mean_response_vscol(sflat_file):
    """
    For an input .fits file, calculates the normalized sigma clipped
    mean flux vs. Col# for a group of Rows returns two arrays for
    the top and bottom section of the CCD
    """
    amc = MaskedCCD(sflat_file)
    ncol = amc.amp_geom.nx
    sensor_type = amc.amp_geom.vendor.lower()
    imaging = amc.amp_geom.imaging
    # use 200 rows close to the amplifier
    row_lo = 10
    row_hi = 210

    # top row
    averow_top = np.zeros(ncol*8)
    for i_amp in range(1, 8+1):
        # Segments 10-17
        anamp = imutils.trim(amc[i_amp], imaging=imaging)
        anamp_im = anamp.getImage()
        anamp_arr = anamp_im.getArray()

        # use a robust mean
        anamp_meanbyrow, _, _ \
            = stats.sigma_clipped_stats(anamp_arr[row_lo:row_hi, :], axis=0)

        # normalize
        nmean_byrow = anamp_meanbyrow/np.median(anamp_meanbyrow)

        lopix = 0 + (i_amp-1)*ncol
        hipix = ncol + (i_amp-1)*ncol
        averow_top[lopix:hipix] = np.flip(nmean_byrow)

    # bot row
    averow_bot = np.zeros((ncol*8))
    for j_amp in range(16, 8, -1):
        # Segments 00-07
        # i_amp goes from 1 to 8, in order of increasing Yccs
        i_amp = 17 - j_amp
        anamp = imutils.trim(amc[j_amp], imaging=imaging)
        anamp_im = anamp.getImage()
        anamp_arr = anamp_im.getArray()

        # use a robust mean
        anamp_meanbyrow, _, _ \
            = stats.sigma_clipped_stats(anamp_arr[row_lo:row_hi, :], axis=0)

        # normalize
        nmean_byrow = anamp_meanbyrow/np.median(anamp_meanbyrow)

        lopix = 0 + (i_amp-1)*ncol
        hipix = ncol + (i_amp-1)*ncol
        if sensor_type == 'e2v':
            averow_bot[lopix:hipix] = nmean_byrow
        elif sensor_type == 'itl':
            averow_bot[lopix:hipix] = np.flip(nmean_byrow)

    # analyze the gaps between amplifiers for Divisidero Tearing, and
    # find the max(abs) deviation in the +-2 columns at the boundaries
    max_divisidero_tearing = []    # 14 entries per CCD
    for k in range(1, 7+1):
        collo = ncol*k - 2  # 2nd to last column in Amplifier
        max_divisidero = np.max(np.abs(averow_top[collo:collo+4] - 1.0))    # +-2 columns
        max_divisidero_tearing.append(max_divisidero)

    for k in range(1, 7+1):
        collo = ncol*k - 2  # 2nd to last column in Amplifier
        max_divisidero = np.max(np.abs(averow_bot[collo:collo+4] - 1.0))    # +-2 columns
        max_divisidero_tearing.append(max_divisidero)

    return averow_top, averow_bot, max_divisidero_tearing


def ana_divisidero_tearing(sflat_files, raft_unit_id, run):
    """
    Analyze a raft of corrected super-flats for Divisidero Tearing.

    Parameters
    ----------
    sflat_files: dict
        Dictionary of lists of single CCD superflat files, keyed by
        slot name.
    raft_unit_id: str
        Raft unit id, e.g., 'LCA-11021_RTM-019'
    run: str
        Run number
    """
    amp_geom = makeAmplifierGeometry(sflat_files['S00'][0])
    ncol = amp_geom.nx

    # make x pixel values
    xpixval = np.arange(ncol*8)

    # dmslotorder
    dmslots = ['S20', 'S21', 'S22', 'S10', 'S11', 'S12', 'S00', 'S01', 'S02']

    # get row averages
    avedict = {}
    for slot in dmslots:
        avedict[slot] = normed_mean_response_vscol(sflat_files[slot][0])

    # make a summary plot
    f = plt.figure(figsize=(20, 20))
    outer = gridspec.GridSpec(3, 3, wspace=0.3, hspace=0.3)

    nskip_edge = 20

    for i, slot in enumerate(dmslots):
        inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[i],
                                                 wspace=0.1, hspace=0.0)
        for j in range(2):

            # use max of max_divisidero_tearing to set the range of plots
            max_divisidero = avedict[slot][2]
            plot_range = np.max(max_divisidero[j*7:j*7+8])

            ax = plt.Subplot(f, inner[j])
            ax.plot(xpixval[nskip_edge:ncol*8 - nskip_edge],
                    avedict[slot][j][nskip_edge:ncol*8
                                     - nskip_edge])
            ax.set_xlabel('Col #')
            ax.set_ylim(1.-plot_range, 1.+plot_range)
            for k in range(1, 8):
                ax.axvline(x=ncol*k, color='red', ls='--', alpha=0.2)
            if j == 0:
                ax.text(0.025, 0.9, '%s' % (slot), transform=ax.transAxes)
                ax.text(0.825, 0.05, 'Seg 10-17', transform=ax.transAxes)
            elif j == 1:
                ax.text(0.825, 0.05, 'Seg 00-07', transform=ax.transAxes)

            f.add_subplot(ax)

    plt.suptitle('Run %s %s' % (str(run), raft_unit_id), fontsize=36)
