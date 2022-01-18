"""
@brief Fit exponential model to mean serial overscans.

@author A. Snyder <snyder18@stanford.edu>
"""
from __future__ import print_function
from __future__ import absolute_import

import numpy as np
import copy
from astropy.io import fits
from scipy.optimize import curve_fit

from lsst.eotest.fitsTools import fitsWriteto, fitsTableFactory
from .MaskedCCD import MaskedCCD
from .AmplifierGeometry import parse_geom_kwd

class OverscanResults(object):

    def __init__(self, all_amps):

        self.output = fits.HDUList()
        self.output.append(fits.PrimaryHDU())
        self.all_amps = all_amps

        self.column_mean = {amp : [] for amp in all_amps}
        self.column_variance = {amp : [] for amp in all_amps}
        self.row_mean = {amp : [] for amp in all_amps}
        self.row_variance = {amp : [] for amp in all_amps}
        self.serial_overscan_noise = {amp : [] for amp in all_amps}
        self.parallel_overscan_noise = {amp : [] for amp in all_amps}
        self.flatfield_signal = {amp : [] for amp in all_amps}
        self.exptime = []   # Frame exposure time
        self.seqnum = []    # Observation sequence number
        self.dayobs = []    # Observation date, e.g., 20211204

    def process_image(self, ccd, gains):
        """Process an image."""
        self.exptime.append(ccd.md.get('EXPTIME'))
        self.seqnum.append(ccd.md.get('SEQNUM'))
        try:
            self.dayobs.append(ccd.md.get('DAYOBS'))
        except KeyError:
            self.dayobs.append(0)
        for amp in self.all_amps:
            image = ccd.bias_subtracted_image(amp)

            datasec = ccd.amp_geom[1]['DATASEC']
            amp_geom = parse_geom_kwd(datasec)
            xmin = amp_geom['xmin']
            xmax = amp_geom['xmax']
            ymin = amp_geom['ymin']
            ymax = amp_geom['ymax']

            imarr = image.getImage().getArray()*gains[amp]

            column_mean = np.mean(imarr[ymin-1:ymax, xmax-1:], axis=0)
            column_variance = np.var(imarr[ymin-1:ymax, xmax-1:], axis=0)
            row_mean = np.mean(imarr[ymax-1:, xmin-1:xmax], axis=1)
            row_variance = np.var(imarr[ymax-1:, xmin-1:xmax], axis=1)
            serial_overscan_noise = np.mean(np.std(imarr[ymin-1:ymax, xmax+2:], axis=1))
            parallel_overscan_noise = np.mean(np.std(imarr[ymax+2:, xmin-1:xmax], axis=1))
            flatfield_signal = np.mean(imarr[ymin-1:ymax, xmin-1:xmax])
                
            self.column_mean[amp].append(column_mean)
            self.column_variance[amp].append(column_variance)
            self.row_mean[amp].append(row_mean)
            self.row_variance[amp].append(row_variance)
            self.serial_overscan_noise[amp].append(serial_overscan_noise)
            self.parallel_overscan_noise[amp].append(parallel_overscan_noise)
            self.flatfield_signal[amp].append(flatfield_signal)

        self.output[0].header['DATASEC'] = datasec

    def build_output_dict(self):
        """Export the results as a dictionary of dictionaries"""
        out_dict = {}
        for amp in self.all_amps:
            extname = 'Amp{0:02d}'.format(amp)
            out_dict[extname] = dict(COLUMN_MEAN=self.column_mean[amp],
                                     COLUMN_VARIANCE=self.column_variance[amp],
                                     ROW_MEAN=self.row_mean[amp],
                                     ROW_VARIANCE=self.row_variance[amp],
                                     FLATFIELD_SIGNAL=self.flatfield_signal[amp],
                                     SERIAL_OVERSCAN_NOISE=self.serial_overscan_noise[amp],
                                     PARALLEL_OVERSCAN_NOISE=self.parallel_overscan_noise[amp])
        return out_dict

    def write_results(self, outfile):
        """Export results as a FITs file."""
        exptime = np.asarray(self.exptime)
        seqnum = np.asarray(self.seqnum)
        dayobs = np.asarray(self.dayobs)

        for amp in self.all_amps:
            extname = 'Amp{0:02d}'.format(amp)
            ncols = len(self.column_mean[amp][0])
            nrows = len(self.row_mean[amp][0])

            idx = np.argsort(self.flatfield_signal[amp])

            cols = [fits.Column('COLUMN_MEAN', format='{0}E'.format(ncols),
                                unit='e-', array=np.asarray(self.column_mean[amp])[idx, :]),
                    fits.Column('COLUMN_VARIANCE', format='{0}E'.format(ncols),
                                unit='e-', array=np.asarray(self.column_variance[amp])[idx, :]),
                    fits.Column('ROW_MEAN', format='{0}E'.format(nrows),
                                unit='e-', array=np.asarray(self.row_mean[amp])[idx, :]),
                    fits.Column('ROW_VARIANCE', format='{0}E'.format(nrows),
                                unit='e-', array=np.asarray(self.row_variance[amp])[idx, :]),
                    fits.Column('FLATFIELD_SIGNAL', format='E', unit='e-',
                                array=np.asarray(self.flatfield_signal[amp])[idx]),
                    fits.Column('SERIAL_OVERSCAN_NOISE', format='E', unit='e-',
                                array=np.asarray(self.serial_overscan_noise[amp])[idx]),
                    fits.Column('PARALLEL_OVERSCAN_NOISE', format='E', unit='e-',
                                array=np.asarray(self.parallel_overscan_noise[amp])[idx]),
                    fits.Column('EXPOSURE', format='E', unit='seconds', array=exptime[idx]),
                    fits.Column('SEQNUM', format='J', unit='None', array=seqnum[idx]),
                    fits.Column('DAYOBS', format='J', unit='None', array=dayobs[idx])
            ]

            self.output.append(fitsTableFactory(cols))
            self.output[-1].name = extname

        self.output[0].header['NAMPS'] = len(self.all_amps)
        fitsWriteto(self.output, outfile, overwrite=True, checksum=True)
