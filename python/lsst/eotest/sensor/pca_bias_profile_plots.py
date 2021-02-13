import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from .ccd_bias_pca import CCD_bias_PCA
from .overscan_frame import make_overscan_frame
from .AmplifierGeometry import makeAmplifierGeometry


__all__ = ['pca_bias_profile_plots']


def plot_imarr(imarr, vmin=-10, vmax=10):
    plt.imshow((imarr - imarr.mean()).T, vmin=vmin, vmax=vmax, origin='lower')
    plt.colorbar()


def pca_bias_profile_plots(raw_file, amp, pca_bias_files, suffix=''):
    amp_geom = makeAmplifierGeometry(raw_file)
    ccd_pcas = CCD_bias_PCA.read_model(*pca_bias_files)
    if not hasattr(ccd_pcas, 'nx'):
        ccd_pcas.nx = 10
        ccd_pcas.ny = 10

    with fits.open(raw_file) as hdus:
        raft = hdus[0].header['RAFTBAY']
        sensor = hdus[0].header['CCDSLOT']
        Run = hdus[0].header['RUNNUM']
        seqnum = hdus[0].header['SEQNUM']
        datasec = hdus[1].header['DATASEC']
        prescan = int(datasec.strip('[]').split(':')[0]) - 1
    bias_file = f'bias_model_{raft}_{sensor}_{Run}_{seqnum:06d}_median.fits'
    residuals_file = f'residuals_{raft}_{sensor}_{Run}_{seqnum:06d}_median.fits'
    ccd_pcas.make_bias_frame(raw_file, bias_file,
                             residuals_file=residuals_file)
    overscan_file = f'overscan_model_{raft}_{sensor}_{Run}_{seqnum:06d}.fits'
    make_overscan_frame(raw_file, outfile=overscan_file)

    with fits.open(raw_file) as raw, fits.open(bias_file) as bias,\
         fits.open(overscan_file) as oscan:
        nrows = 6
        row_height = 3
        alpha = 0.5
        print(amp)
        title = f'Run {Run}, {raft}_{sensor}, SEQNUM {seqnum}, amp {amp}'
        fig = plt.figure(figsize=(10, nrows*row_height))
        fig.add_subplot(nrows, 1, 1)
        plot_imarr(raw[amp].data)
        plt.title('raw bias image')
        fig.add_subplot(nrows, 1, 2)
        plot_imarr(bias[amp].data)
        plt.title('PCA bias model')
        fig.add_subplot(nrows, 1, 3)
        plot_imarr(raw[amp].data - bias[amp].data)
        plt.title('raw - bias model')
        fig.add_subplot(nrows, 1, 4)
        plot_imarr(oscan[amp].data)
        plt.title('overscan-based image')
        fig.add_subplot(nrows, 1, 5)
        plot_imarr(raw[amp].data - oscan[amp].data)
        plt.title('raw - overscan-based image')
        plt.suptitle('\n'.join((title, os.path.basename(raw_file))))


        # Plot raw, bias model, and overscan profiles in the serial
        # direction using the pixels in the parallel overscan region.
        rows = slice(amp_geom.parallel_overscan.getMinY(),
                     amp_geom.parallel_overscan.getMaxY())
        fig.add_subplot(nrows, 2, nrows*2-1)
        oscan_profile = np.mean(oscan[amp].data[rows, :], axis=0)[prescan:]
        plt.plot(oscan_profile, label='overscan region', alpha=alpha)
        plt.plot(np.mean(raw[amp].data, axis=0)[prescan:],
                 label='raw data (full segment)', alpha=alpha)
        bias_profile = np.mean(bias[amp].data[rows, :], axis=0)[prescan:]
        plt.plot(bias_profile, label='bias model (overscan region)',
                 alpha=alpha)
        ymin, ymax = np.min(bias_profile), np.max(bias_profile)
        dy = (ymax - ymin)/5.
        plt.ylim(ymin - dy, ymax + dy)
        plt.legend(fontsize='x-small')
        plt.title('serial direction profiles')

        # Plot raw, bias model, and overscan profiles in the parallel
        # direction using the pixels in the serial overscan region.
        columns = slice(amp_geom.serial_overscan.getMinX(),
                        amp_geom.serial_overscan.getMaxX())
        fig.add_subplot(nrows, 2, nrows*2)
        oscan_profile = np.mean(oscan[amp].data[:, columns], axis=1)
        plt.plot(oscan_profile, label='overscan region',
                 alpha=alpha)
        plt.plot(np.mean(raw[amp].data, axis=1),
                 label='raw data (full segment)', alpha=alpha)
        plt.plot(np.mean(bias[amp].data[:, columns], axis=1),
                 label='bias model (overscan region)', alpha=alpha)
        ymin, ymax = np.min(oscan_profile), np.max(oscan_profile)
        dy = (ymax - ymin)/5.
        plt.ylim(ymin - dy, ymax + dy)
        plt.title('parallel direction profiles')
        plt.savefig(f'{Run}_{raft}_{sensor}_{seqnum}_{amp}_{suffix}_'
                    'bias_model_profiles.png')

        fig1 = plt.figure(figsize=(10, 8))
        plt.hist((raw[amp].data - bias[amp].data)[:, prescan:].ravel(),
                 bins=100, range=(-50, 50), histtype='step',
                 label='raw - bias model')
        plt.hist((raw[amp].data - oscan[amp].data)[:, prescan:].ravel(),
                 bins=100, range=(-50, 50), histtype='step',
                 label='raw - overscan-based image')
        plt.title(title)
        plt.legend(fontsize='x-small')
        plt.xlabel('residuals (ADU/pixel)')
        plt.yscale('log')
        plt.savefig(f'{Run}_{raft}_{sensor}_{seqnum}_{amp}_{suffix}_'
                    'bias_model_residuals_hist.png')
