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

    def __init__(self, outfile=None):

        self.outfile = outfile
        if outfile is None:
            self.output = fits.HDUList()
            self.output.append(fits.PrimaryHDU())
        else:
            self.output = fits.open(self.outfile)

        self.meanrow = {amp : [] for amp in range(1, 17)}
        self.var = {amp : [] for amp in range(1, 17)}
        self.oscan_noise = {amp : [] for amp in range(1, 17)}
        self.flux = {amp : [] for amp in range(1, 17)}
        self.flux_std = {amp : [] for amp in range(1, 17)}

    def process_image(self, ccd, gains):
        """Process an image."""
      
        for amp in range(1, 17):
            image = ccd.bias_subtracted_image(amp)

            datasec = ccd.amp_geom[1]['DATASEC']
            amp_geom = parse_geom_kwd(datasec)
            xmin = amp_geom['xmin']
            xmax = amp_geom['xmax']
            ymin = amp_geom['ymin']
            ymax = amp_geom['ymax']

            imarr = image.getImage().getArray()*gains[amp]

            meanrow = np.mean(imarr[ymin-1:ymax, :], axis=0)
            var = np.var(imarr[ymin-1:ymax, :], axis=0)
            oscan_noise = np.mean(np.std(imarr[ymin-1:ymax, xmax+2:xmax+28], axis=1))
            flux = np.mean(imarr[ymin-1:ymax, xmin-1:xmax])
            flux_std = np.std(imarr[ymin-1:ymax, xmin-1:xmax])
                
            self.meanrow[amp].append(meanrow)
            self.var[amp].append(var)
            self.flux[amp].append(flux)
            self.flux_std[amp].append(flux_std)
            self.oscan_noise[amp].append(oscan_noise)

        self.output[0].header['DATASEC'] = datasec

    def build_output_dict(self):
        """Export the results as a dictionary of dictionaries"""
        out_dict = {}
        for amp in range(1, 17):
            extname = 'Amp{0:02d}'.format(amp)
            out_dict[extname] = dict(MEANROW=self.meanrow[amp],
                                     VAR=self.var[amp],
                                     FLUX=self.flux[amp],
                                     FLUX_STD=self.flux_std[amp],
                                     NOISE=self.noise[amp])
        return out_dict

    def write_results(self, outfile):
        """Export results as a FITs file."""
        
        for amp in range(1, 17):
            extname = 'Amp{0:02d}'.format(amp)
            nrows1 = len(self.flux[amp])
            ncols = len(self.meanrow[amp][0])

            meanrow_col = fits.Column('MEANROW', format='{0}E'.format(ncols), unit='e-', 
                                      array=self.meanrow[amp])
            var_col = fits.Column('VAR', format='{0}E'.format(ncols), unit='e-',
                                  array=self.var[amp])
            flux_col = fits.Column('FLUX', format='E', unit='e-', array=self.flux[amp])
            flux_std_col = fits.Column('FLUX_STD', format='E', unit='e-', 
                                       array=self.flux_std[amp])
            noise_col = fits.Column('NOISE', format='E', unit='e-', array=self.noise[amp])

            try:
                #
                # Append new rows if HDU for this segment already exists
                #
                table_hdu = self.output[extname]
                row0 = table_hdu.header['NAXIS2']
                nrows = row0+nrows1
                table_hdu = fitsTableFactory(table_hdu.data, nrows=nrows)
                table_hdu.data['MEANROW'][row0:] = meanrow_col
                table_hdu.data['VAR'][row0:] = var_col
                table_hdu.data['FLUX'][row0:] = flux_col
                table_hdu.data['FLUX_STD'][row0:] = flux_std_col
                table_hdu.data['NOISE'][row0:] = noise_col
                table_hdu.name = extname
                self.output[extname] = table_hdu
            except KeyError:
                self.output.append(fitsTableFactory([meanrow_col, 
                                                     var_col,
                                                     flux_col,
                                                     flux_std_col,
                                                     noise_col]))
                self.output[-1].name = extname

        self.output[0].header['NAMPS'] = 16
        fitsWriteto(self.output, outfile, overwrite=True, checksum=True)
