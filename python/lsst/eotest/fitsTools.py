import warnings
import numpy as np
import astropy.io.fits as fits


def fitsTableFactory(*args, **kwds):
    """Wrap the binary table hdu generation code."""
    my_table = fits.BinTableHDU.from_columns(*args, **kwds)
    return my_table


def fitsWriteto(hdulist, outfile, **kwds):
    """Silence warnings about over-writing a file."""
    warnings.resetwarnings()
    warnings.filterwarnings('ignore', category=UserWarning, append=True)
    hdulist.writeto(outfile, **kwds)
    warnings.resetwarnings()
    warnings.filterwarnings('always', category=UserWarning, append=True)


if __name__ == '__main__':
    fits_cols = [fits.Column(name='foo', format='I', unit='None',
                             array=np.zeros(10))]
    foo = fitsTableFactory(fits_cols)
