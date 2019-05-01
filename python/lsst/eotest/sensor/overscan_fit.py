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

def make_joint_model(flux, N):
    """Make high flux exponential plus CTI model."""
    
    def joint_model(x, A, tau, cti):
        
        result = A*np.exp(-x/tau) + (cti**x)*N*flux
        
        return result
        
    return joint_model

class OverscanFit(object):

    def __init__(self, num_oscan_pixels=10, minflux=30000., maxflux=150000., 
                 outfile=None):

        self.num_oscan_pixels = num_oscan_pixels
        if maxflux <= minflux:
            raise ValueError("maxflux must be greater than minflux")
        self.minflux = minflux
        self.maxflux = maxflux
        self.outfile = outfile
        if outfile is None:
            self.output = fits.HDUList()
            self.output.append(fits.PrimaryHDU())
        else:
            self.output = fits.open(self.outfile)

        self.meanrow = {amp : [] for amp in range(1, 17)}
        self.noise = {amp : [] for amp in range(1, 17)}
        self.flux = {amp : [] for amp in range(1, 17)}
        self.flux_std = {amp : [] for amp in range(1, 17)}
        self.signal = {amp : [] for amp in range(1, 17)}
        self.signal_std = {amp : [] for amp in range(1, 17)}
        self.tau = {amp : [] for amp in range(1, 17)}
        self.tau_std = {amp : [] for amp in range(1, 17)}
        self.cti = {amp : [] for amp in range(1, 17)}
        self.cti_std = {amp : [] for amp in range(1, 17)}

    def process_image(self, ccd, gains):
        """Process an image."""

        for amp in range(1, 17):
            image = ccd.bias_subtracted_image(amp)

            ## get xmin, xmax, ymin, ymax here
            amp_geom = parse_geom_kwd(ccd.amp_geom[1]['DATASEC'])
            xmin = amp_geom['xmin']
            xmax = amp_geom['xmax']
            ymin = amp_geom['ymin']
            ymax = amp_geom['ymax']

            imarr = image.getImage().getArray()*gains[amp]

            meanrow = np.mean(imarr[ymin-1:ymax, :], axis=0)
            noise = np.mean(np.std(imarr[ymin-1:ymax, xmax+2:xmax+28], axis=1))
            flux = np.mean(imarr[ymin-1:ymax, xmin-1:xmax])
            flux_std = np.std(imarr[ymin-1:ymax, xmin-1:xmax])
            signal = np.nan
            tau = np.nan
            cti = np.nan
            signal_std = np.nan
            tau_std = np.nan
            cti_std = np.nan
        
            ## Perform overscan exponential function fit
            if self.minflux <= flux <= self.maxflux:
                y = copy.deepcopy(meanrow[xmax:xmax+self.num_oscan_pixels])
                x = np.arange(1, y.shape[0]+1)

                try:
                    fit_params, fit_covar = curve_fit(make_joint_model(flux, xmax),
                                                      x, y, p0=(50, 1.25, 1E-6))
                except RuntimeError:
                    print("Error for fit: flux = {0:.1f}, Skipping...".format(flux))
                else:
                    signal = fit_params[0]
                    tau = fit_params[1]
                    cti = fit_params[2]
                    signal_std = np.sqrt(fit_covar[0, 0])
                    tau_std = np.sqrt(fit_covar[1, 1])
                    cti_std = np.sqrt(fit_covar[2, 2])
                
            self.meanrow[amp].append(meanrow)
            self.flux[amp].append(flux)
            self.flux_std[amp].append(flux_std)
            self.noise[amp].append(noise)
            self.signal[amp].append(signal)
            self.tau[amp].append(tau)
            self.cti[amp].append(cti)
            self.signal_std[amp].append(signal_std)
            self.tau_std[amp].append(tau_std)
            self.cti_std[amp].append(cti_std)

    def write_results(self, outfile):
        """Export results as a FITs file."""
        
        for amp in range(1, 17):
            extname = 'Amp{0:02d}'.format(amp)
            nrows1 = len(self.flux[amp])
            ncols = len(self.meanrow[amp][0])

            meanrow_col = fits.Column('MEANROW', format='{0}E'.format(ncols), unit='e-', 
                                      array=self.meanrow[amp])
            flux_col = fits.Column('FLUX', format='E', unit='e-', array=self.flux[amp])
            flux_std_col = fits.Column('FLUX_STD', format='E', unit='e-', 
                                       array=self.flux_std[amp])
            noise_col = fits.Column('NOISE', format='E', unit='e-', array=self.noise[amp])
            signal_col = fits.Column('SIGNAL', format='E', unit='e-', array=self.signal[amp])
            signal_std_col = fits.Column('SIGNAL_STD', format='E', unit='e-', 
                                         array=self.signal_std[amp])
            tau_col = fits.Column('TAU', format='E', unit='None', array=self.tau[amp])
            tau_std_col = fits.Column('TAU_STD', format='E', unit='None', 
                                      array=self.tau_std[amp])
            cti_col = fits.Column('CTI', format='E', unit='None', array=self.cti[amp])
            cti_std_col = fits.Column('CTI_STD', format='E', unit='None', 
                                      array=self.cti_std[amp])

            try:
                #
                # Append new rows if HDU for this segment already exists
                #
                table_hdu = self.output[extname]
                row0 = table_hdu.header['NAXIS2']
                nrows = row0+nrows1
                table_hdu = fitsTableFactory(table_hdu.data, nrows=nrows)
                table_hdu.data['MEANROW'][row0:] = meanrow_col
                table_hdu.data['FLUX'][row0:] = flux_col
                table_hdu.data['FLUX_STD'][row0:] = flux_std_col
                table_hdu.data['NOISE'][row0:] = noise_col
                table_hdu.data['SIGNAL'][row0:] = signal_col
                table_hdu.data['SIGNAL_STD'][row0:] = signal_std_col
                table_hdu.data['TAU'][row0:] = tau_col
                table_hdu.data['TAU_STD'][row0:] = tau_std_col
                table_hdu.data['CTI'][row0:] = cti_col
                table_hdu.data['CTI_STD'][row0:] = cti_std_col
                table_hdu.name = extname
                self.output[extname] = table_hdu
            except KeyError:
                self.output.append(fitsTableFactory([meanrow_col, 
                                                     flux_col,
                                                     flux_std_col,
                                                     noise_col,
                                                     signal_col,
                                                     signal_std_col,
                                                     tau_col,
                                                     tau_std_col,
                                                     cti_col,
                                                     cti_std_col]))
                self.output[-1].name = extname

        self.output[0].header['NAMPS'] = 16
        fitsWriteto(self.output, outfile, overwrite=True, checksum=True)
