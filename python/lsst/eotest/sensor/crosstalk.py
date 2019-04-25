from __future__ import print_function
from __future__ import absolute_import
import numpy as np
from copy import deepcopy
from astropy.io import fits
from scipy.ndimage.filters import gaussian_filter

from lsst.eotest.fitsTools import fitsWriteto

def is_valid_aggressor(ccd, amp, ay, ax, threshold=100000., r=50):
    """Evaluate if identified spot is valid aggressor."""    
    print(ay, ax)
    imarr = ccd[amp].getImage().getArray()
    ny, nx = imarr.shape
    y, x = np.ogrid[-ay:ny-ay, -ax:nx-ax]
    mask = x*x + y*y >= r*r
    spot_arr = np.ma.MaskedArray(imarr, mask)
    print(np.mean(spot_arr))

    return np.mean(spot_arr) > threshold
    
def find_aggressors(ccd, threshold=100000, gf_sigma=50):
    """Determine amplifier location of aggressor spots."""
    candidate_list = []
    for amp in ccd.keys():
        blurred = gaussian_filter(ccd[amp].getImage().getArray(), gf_sigma)
        y, x = np.unravel_index(blurred.argmax(), blurred.shape)
        
        if is_valid_aggressor(ccd, amp, y, x, threshold=threshold):
            candidate_list.append((amp, y, x))

    return candidate_list

def make_stamp(ccd, amp, y, x, l=200):
    """Get postage stamp for crosstalk calculations."""

    imarr = ccd.unbiased_and_trimmed_image(amp).getImage().getArray()
    maxy, maxx = imarr.shape

    ## Truncate at amplifier edges
    if y-l//2 < 0: 
        y0 = 0
    else:
        y0 = y-l//2
    if y+l//2 >= maxy: 
        y1 = maxy
    else:
        y1 = y+l//2
    if x-l//2 < 0: 
        x0 = 0
    else:
        x0 = x-l//2
    if x+l//2 >= maxx: 
        x1 = maxx
    else:
        x1 = x+l//2

    return deepcopy(imarr[y0:y1, x0:x1])

def crosstalk_model(coefficients, aggressor_array):
    """Create a crosstalk victim postage stamp from model."""

    ny, nx = aggressor_array.shape
    crosstalk_signal = coefficients[0]
    bias = coefficients[1]
    tilty = coefficients[2]
    tiltx = coefficients[3]

    Y, X = np.mgrid[:ny, :nx]
    model = crosstalk_signal*aggressor_array + tilty*Y + tiltx*X + bias
    return model

def crosstalk_model_fit(aggressor_stamp, victim_stamp, num_iter=5, nsig=2.0):
    """Perform a crosstalk model fit for given  aggressor and victim stamps."""

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

    return np.append(coefficients[0], np.sqrt(covar.diagonal()))
        
class CrosstalkMatrix():

    def __init__(self, aggressor_sensor_id, victim_sensor_id=None, filename=None, namps=16):
        self.header = fits.Header()
        self.header.set('SENSOR1', aggressor_sensor_id)
        if victim_sensor_id is not None:
            self.header.set('SENSOR2', victim_sensor_id)
        else:
            self.header.set('SENSOR2', aggressor_sensor_id)
        self.filename = filename
        self.namps = namps
        self._set_matrix()
        if self.filename is not None:
            self._read_matrix()

    def set_row(self, aggressor_amp, row):
        """Set matrix row from results dictionary"""
        for victim_amp in row.keys():
            self.matrix[:, aggressor_amp-1, victim_amp-1] = row[victim_amp]

    def _set_matrix(self):
        """Initialize crosstalk matrix as NaNs."""
        self.matrix = np.zeros((8, self.namps, self.namps), dtype=np.float)
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
        for i in range(8):
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
        xtalk_hdu = fits.PrimaryHDU(self.matrix[0,:,:], header=self.header)
        bias_hdu = fits.ImageHDU(self.matrix[1,:,:], name='BIAS')
        tilty_hdu = fits.ImageHDU(self.matrix[2,:,:], name='TILT_Y')
        tiltx_hdu = fits.ImageHDU(self.matrix[3,:,:], name='TILT_X')
        xtalkerr_hdu = fits.ImageHDU(self.matrix[4,:,:], name='SIGMA_XTALK')
        biaserr_hdu = fits.ImageHDU(self.matrix[5,:,:], name='SIGMA_BIAS')
        tiltyerr_hdu = fits.ImageHDU(self.matrix[6,:,:], name='SIGMA_TILT_Y')
        tiltxerr_hdu = fits.ImageHDU(self.matrix[7,:,:], name='SIGMA_TILT_X')
        
        output = fits.HDUList([xtalk_hdu, bias_hdu, tilty_hdu, tiltx_hdu, 
                               xtalkerr_hdu, biaserr_hdu, tiltyerr_hdu, tiltxerr_hdu])
        fitsWriteto(output, outfile, overwrite=overwrite)

    def write_yaml(self, outfile, overwrite=True):
        """Write crosstalk results to a YAML file."""
        raise NotImplementedError



        
    
    
                
            

        
            
            
        
        

            
            
        
        
