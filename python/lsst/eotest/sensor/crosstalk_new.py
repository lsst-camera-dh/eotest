from __future__ import print_function
from __future__ import absolute_import
import numpy as np
from astropy.io import fits
from scipy.ndimage.filters import gaussian_filter

from lsst.eotest.fitsTools import fitsWriteto

def is_valid_aggressor(ccd, amp, ay, ax, threshold=100000, r=50):
    """Determine if footprint is a valid spot"""    
    imarr = ccd[amp].getImage().getArray()
    ny, nx = imarr.shape
    y, x = np.ogrid[-ay:ny-ay, -ax:nx-ax]
    mask = x*x + y*y >= r*r
    spot_arr = np.ma.MaskedArray(imarr, mask)

    return np.mean(spot_array) > threshold
    
def aggressors(ccd, threshold=100000, gf_sigma=50):
    """Locate aggressor amplifiers"""
    candidate_list = []
    for amp in ccd.keys():
        blurred = gaussian_filter(ccd[amp].getImage().getArray(), gf_sigma)
        y, x = np.unravel_index(blurred.argmax(), blurred.shape)
        
        if is_valid_aggressor(ccd, amp, y, x, threshold=threshold)
            candidate_list.append((amp, y, x))

    return candidate_list

def crosstalk_model(coefficients, aggressor_array):
    """Make a model of a crosstalk victim postage stamp."""

    ny, nx = aggressor_array.shape
    crosstalk_signal = coefficients[0]
    bias = coefficients[1]
    tilty = coefficients[2]
    tiltx = coefficients[3]

    Y, X = np.mgrid[:ny, :nx]
    model = crosstalk_signal*aggressor_array + tilty*Y + tiltx*X + bias
    return model

def crosstalk_model_fit(aggressor_stamp, victim_stamp, num_iter=10, nsig=3.0):

    coefficients = np.asarray([[0,0,0,0]])
    victim_array = np.ma.masked_invalid(victim_stamp)
    mask = np.ma.getmask(victim_array)

    for i in range(num_iter):
        #
        # Mask outliers using residual
        #
        model = np.ma.masked_where(mask, crosstalk_model(coefficients[0],
                                                         aggressor_stamp))
        residual = victim_array - model
        res_mean = residual.mean()
        res_std = residual.std()
        victim_array = np.ma.masked_where(np.abs(residual-res_mean) \
                                              > nsig*res_std, victim_stamp)
        mask = np.ma.getmask(victim_array)
        #
        # Construct masked, compressed basis arrays
        #
        ay, ax = aggressor_stamp.shape
        bias = np.ma.masked_where(mask, np.ones((ay, ax))).compressed()
        Y, X = np.mgrid[:ay, :ax]
        Y = np.ma.masked_where(mask, Y).compressed()
        X = np.ma.masked_where(mask, X).compressed()
        aggressor_array = np.ma.masked_where(mask, aggressor_stamp).compressed()
        #
        # Perform least-squares minimization
        #
        b = victim_array.compressed()/victim_array.std()
        A = np.vstack([aggressor_array, bias, Y, X]).T/victim_array.std()
        coefficients = np.linalg.lstsq(A, b)
        covar = np.matrix(np.dot(A.T, A)).I

    return np.append(crosstalk_results[0], np.sqrt(covar[0, 0]))

def calculate_crosstalk(aggressor_ccd, victim_ccd, threshold=100000,
                        num_iter=10, nsig=3.0, crosstalk_matrix=None):
    """Calculate crosstalk values from a single crosstalk image"""
    #
    # Initialize new crosstalk matrix if none given
    #
    if crosstalk_matrix is None:
        crosstalk_matrix = CrosstalkMatrix()
    #
    # Find aggressor amplifiers
    #
    aggressor_list = aggressors(aggressor_ccd, threshold)
    for aggressor_amp, y, x in aggressor_amp_list:

        aggressor_stamp = stamp(aggressor_ccd, y, x)
        row = dict()
        for amp in victim_ccd.keys():
            victim_stamp = stamp(victim_ccd, y, x)
            try:
                row[amp] = crosstalk_model_fit(aggressor_stamp, victim_stamp,
                                               num_iter=num_iter, nsig=nsig)
            except: # check what error is thrown for failed lstsq 
                print("Error extracting crosstalk signal.")
                print("Skipping")
        crosstalk_matrix.set_row(aggressor_amp, row)

    return crosstalk_matrix

def make_crosstalk_matrix(aggressor_infiles, victim_infiles=None, 
                          treshold=100000, nsig=3.0):
    #
    # May change if wish not to save to memory
    #
    crosstalk_matrix = CrosstalkMatrix()
    for n, aggressor_infile in enumerate(aggressor_infiles):
        aggressor_ccd = MaskedCCD(aggressor_infile)
        #
        # Set victim ccd if provided
        #
        if victim_infiles is None:
            victim_ccd = MaskedCCD(aggressor_infile)
        else:
            try:
                victim_ccd = MaskedCCD(victim_infiles[n])
            except IndexError:
                print("Missing victim image file.")
                print("Skipping")
                continue

        crosstalk_matrix = calculate_crosstalk(aggressor_ccd, victim_ccd, 
                                               threshold=threshold, 
                                               nsig=nsig, crosstalk_matrix)

    return crosstalk_matrix
        
class CrosstalkMatrix():

    def __init__(self, filename=None, namps=16):
        self.filename = filename
        self.namps = namps
        self._set_matrix()
        if self.filename is not None:
            self._read_matrix()

    def set_row(self, agg, row):
        """Set matrix row from results dictionary"""
        for vic in row.keys():
            self.matrix[:, agg, vic] = row[vic]

    def _set_matrix(self):
        """Initialize crosstalk matrix as NaNs."""
        self.matrix = np.zeros((5, self.namps, self.namps), dtype=np.float)
        self.matrix[:] = np.nan

    def _read_matrix(self):
        """Read crosstalk matrix from file."""
        if self.filename[-5:] == '.fits':
            self._read_fits_matrix()
        elif self.filename[-5:] == '.yaml':
            self._read_yaml_matrix()
        else:
            raise ValueError('Crosstalk matrix file must be FITS or YAML filetype')

    def _read_fits_matrix(self):
        """Read crosstalk results from FITS file."""        
        hdulist = fits.open(self.filename)
        for i in range(5):
            self.matrix[i,:,:] = hdulist[i].data

    def _read_yaml_matrix(self):
        """Read crosstalk results from a YAML file."""
        raise NotImplementedError

    def write_fits(self, outfile=None, overwrite=True):
        """Write crosstalk results to FITS file."""
        if outfile is None:
            outfile = self.filename
        else:
            self.filename = outfile
        #
        # Save matrix results into separate HDUs
        #
        xtalk_hdu = fits.PrimaryHDU(self.matrix[0,:,:])
        bias_hdu = fits.PrimaryHDU(self.matrix[1,:,:], name='BIAS')
        tilty_hdu = fits.ImageHDU(self.matrix[2,:,:], name='TILT_Y')
        tiltx_hdu = fits.ImageHDU(self.matrix[3,:,:], name='TILT_X')
        error_hdu = fits.ImageHDU(self.matrix[4,:,:], name='SIGMA_XTALK')
        
        output = fits.HDUList([xtalk_hdu, bias_hdu, tilty_hdu, tiltx_hdu, error_hdu])
        fitsWriteto(output, outfile, overwrite=overwrite)

    def write_yaml(self, outfile, overwrite=True):
        """Write crosstalk results to a YAML file."""
        raise NotImplementedError

    def _add_(self, other):
        """Addition of two matrices (error adds in quadrature)"""
        result = CrosstalkMatrix()
        result.matrix[:4,:,:] = self.matrix[:4,:,:] + other.matrix[:4,:,:]
        result.matrix[4,:,:] = np.sqrt(np.square(self.matrix[4,:,:]) \
                                           + np.square(other.matrix[4,:,:]))
        return result

    def _sub_(self, other):
        """Subtraction of two matrices (error adds in quadrature)"""
        result = CrosstalkMatrix()
        result.matrix[:4,:,:] = self.matrix[:4,:,:] - other.matrix[:4,:,:]
        result.matrix[4,:,:] = np.sqrt(np.square(self.matrix[4,:,:]) \
                                           + np.square(other.matrix[4,:,:]))
        return result



        
    
    
                
            

        
            
            
        
        

            
            
        
        
