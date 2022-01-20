"""
Code to perform Divisadero tearing analysis.  This is slightly
revised code originally from Aaron Roodman. See LSSTTD-1440.
"""
import numpy as np
import numpy.ma as ma
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from astropy import stats
import lsst.eotest.image_utils as imutils
import lsst.eotest.sensor as sensorTest


__all__ = ['ana_divisidero_tearing']


def normed_mean_response_vscol(sflat_file, mask_files=()):
    """
    For an input .fits file, calculates the normalized sigma clipped
    mean flux vs. Col# for a group of Rows returns two arrays for
    the top and bottom section of the CCD
    """
    amc = sensorTest.MaskedCCD(sflat_file, mask_files=mask_files)
    amps = imutils.allAmps(sflat_file)
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
        anamp_arr = ma.masked_array(anamp_im.getArray(),
                                    mask=anamp.getMask().array)

        # use a robust mean
        anamp_meanbyrow, _, _ \
            = stats.sigma_clipped_stats(anamp_arr[row_lo:row_hi, :], axis=0)

        # normalize
        nmean_byrow = anamp_meanbyrow/np.median(anamp_meanbyrow)

        # fit nmean_byrow to a line and divide that line out
        nedge = 25
        x = np.arange(nmean_byrow.shape[0])
        y = nmean_byrow
        cpoly = np.polyfit(x[nedge:-nedge],y[nedge:-nedge],deg=1)
        yfit = cpoly[1] + cpoly[0]*x
        nmean_byrow = y/yfit

        lopix = 0 + (i_amp-1)*ncol
        hipix = ncol + (i_amp-1)*ncol
        averow_top[lopix:hipix] = np.flip(nmean_byrow)

    # bot row
    averow_bot = np.zeros((ncol*8))
    for j_amp in range(16, 8, -1):
        if j_amp not in amps:
            continue
        # Segments 00-07
        # i_amp goes from 1 to 8, in order of increasing Yccs
        i_amp = 17 - j_amp
        anamp = imutils.trim(amc[j_amp], imaging=imaging)
        anamp_im = anamp.getImage()
        anamp_arr = ma.masked_array(anamp_im.getArray(),
                                    mask=anamp.getMask().array)

        # use a robust mean
        anamp_meanbyrow, _, _ \
            = stats.sigma_clipped_stats(anamp_arr[row_lo:row_hi, :], axis=0)

        # normalize
        nmean_byrow = anamp_meanbyrow/np.median(anamp_meanbyrow)

        # fit nmean_byrow to a line and divide that line out
        nedge = 25
        x = np.arange(nmean_byrow.shape[0])
        y = nmean_byrow
        cpoly = np.polyfit(x[nedge:-nedge],y[nedge:-nedge],deg=1)
        yfit = cpoly[1] + cpoly[0]*x
        nmean_byrow = y/yfit

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
        max_divisidero = np.nanmax(np.abs(averow_top[collo:collo+4] - 1.0))    # +-2 columns
        if np.isnan(max_divisidero):
            max_divisidero = 0
        max_divisidero_tearing.append(max_divisidero)

    for k in range(1, 7+1):
        if k + 8 not in amps:
            continue
        collo = ncol*k - 2  # 2nd to last column in Amplifier
        max_divisidero = np.nanmax(np.abs(averow_bot[collo:collo+4] - 1.0))    # +-2 columns
        if np.isnan(max_divisidero):
            max_divisidero = 0
        max_divisidero_tearing.append(max_divisidero)

    return averow_top, averow_bot, max_divisidero_tearing


def ana_divisidero_tearing(sflat_files, mask_files, title=None):
    """
    Analyze a raft of corrected super-flats for Divisidero Tearing.

    Parameters
    ----------
    sflat_files: dict
        Dictionary of single CCD superflat files, keyed by slot name.
    mask_files: dict of lists
        Mask files for each CCD, keyed by slot name.
    title: str [None]
        Plot title.
    """
    my_slot = list(sflat_files)[0]
    amp_geom = sensorTest.makeAmplifierGeometry(sflat_files[my_slot])
    ncol = amp_geom.nx

    # make x pixel values
    xpixval = np.arange(ncol*8)

    # dmslotorder
    dmslots = ['S20', 'S21', 'S22', 'S10', 'S11', 'S12', 'S00', 'S01', 'S02']
    if 'SW0' in sflat_files:
        dmslots = 'SW0 SW1 SG0 SG1'.split()

    # get row averages
    avedict = {}
    for slot in dmslots:
        try:
            avedict[slot] = normed_mean_response_vscol(sflat_files[slot],
                                                       mask_files=mask_files[slot])
        except KeyError:
            # This will occur if data from `slot` is not available.
            pass

    # make a summary plot
    f = plt.figure(figsize=(20, 20))
    outer = gridspec.GridSpec(3, 3, wspace=0.3, hspace=0.3)

    nskip_edge = 20

    for i, slot in enumerate(avedict):
        max_divisidero = avedict[slot][2]
        have_wf_sensor = (len(max_divisidero) == 7)
        inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[i],
                                                 wspace=0.1, hspace=0.0)
        for j in range(2):
            if have_wf_sensor and j==1:
                continue
            # use max of max_divisidero_tearing to set the range of plots
            plot_range = np.nanmax(max_divisidero[j*7:j*7+8])

            ax = plt.Subplot(f, inner[j])
            ax.plot(xpixval[nskip_edge:ncol*8 - nskip_edge],
                    avedict[slot][j][nskip_edge:ncol*8 - nskip_edge])
            ax.set_xlabel('Col #')
            try:
                ax.set_ylim(1.-plot_range, 1.+plot_range)
            except ValueError as eobj:
                # plot_range is probably inf or NaN because of bad pixel
                # data for this sensor, so just skip this plot.
                print('ValueError:', str(eobj))
                continue
            for k in range(1, 8):
                ax.axvline(x=ncol*k, color='red', ls='--', alpha=0.2)
            if j == 0 and not have_wf_sensor:
                ax.text(0.025, 0.9, '%s' % (slot), transform=ax.transAxes)
                ax.text(0.825, 0.05, 'Seg 10-17', transform=ax.transAxes)
            elif j == 1 or have_wf_sensor:
                ax.text(0.825, 0.05, 'Seg 00-07', transform=ax.transAxes)

            f.add_subplot(ax)

    plt.suptitle(title, fontsize=36)
    return {slot: avedict[slot][2] for slot in avedict}
