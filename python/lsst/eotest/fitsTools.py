import warnings
import numpy as np
import astropy
from astropy.io import fits


def fitsTableFactory(*args, **kwds):
    """Wrap the binary table hdu generation code."""
    my_table = fits.BinTableHDU.from_columns(*args, **kwds)
    return my_table


def fitsWriteto(hdulist, outfile, compress_images=True, **kwds):
    """
    Silence warnings about over-writing a file and optionally perform
    compression on each ImageHDU.
    """
    # Handle change from clobber to overwrite:
    if 'overwrite' in kwds and astropy.__version__.startswith('1.'):
        kwds['clobber'] = kwds['overwrite']
        del kwds['overwrite']
    warnings.resetwarnings()
    warnings.filterwarnings('ignore', category=UserWarning, append=True)
    warnings.filterwarnings('ignore', category=fits.verify.VerifyWarning,
                            append=True)
    if isinstance(hdulist, fits.HDUList) and compress_images:
        for i in range(len(hdulist)):
            if isinstance(hdulist[i], fits.ImageHDU):
                hdulist[i] = fits.CompImageHDU(data=hdulist[i].data,
                                               header=hdulist[i].header)
    hdulist.writeto(outfile, **kwds)


if __name__ == '__main__':
    fits_cols = [fits.Column(name='foo', format='I', unit='None',
                             array=np.zeros(10))]
    foo = fitsTableFactory(fits_cols)
