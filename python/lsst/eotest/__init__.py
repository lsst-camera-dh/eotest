import os
import lsst.utils
from .version import *

def getVersion():
    if __version__ != "unknown":
        return __version__
    import lsst.eotest
    for item in lsst.eotest.__file__.split(os.path.sep):
        if item.startswith('eotest-'):
            return item[len('eotest-'):]
