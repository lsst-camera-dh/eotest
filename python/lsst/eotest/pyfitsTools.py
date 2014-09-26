import warnings
import numpy as np
import pyfits

def pyfitsTableFactory(*args, **kwds):
    """Silence deprecation warning about pyfits.new_table in pyfits 3.3"""
    warnings.resetwarnings()
    warnings.filterwarnings('ignore', category=UserWarning, append=True)
    my_table = pyfits.new_table(*args, **kwds)
    warnings.resetwarnings()
    warnings.filterwarnings('always', category=UserWarning, append=True)
    return my_table

def pyfitsWriteto(hdulist, outfile, **kwds):
    """Silence warnings about over-writing a file."""
    warnings.resetwarnings()
    warnings.filterwarnings('ignore', category=UserWarning, append=True)
    hdulist.writeto(outfile, **kwds)
    warnings.resetwarnings()
    warnings.filterwarnings('always', category=UserWarning, append=True)

if __name__ == '__main__':
    fits_cols = [pyfits.Column(name='foo', format='I', unit='None',
                               array=np.zeros(10))]
    foo = pyfits_table_factory(fits_cols)
