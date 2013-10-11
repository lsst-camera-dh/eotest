"""
@brief Provide a standard interface for accessing segment gains, either
from user input or by querying the Sensor Test database.

@author J. Chiang <jchiang@.slac.stanford.edu>
"""
import os
from image_utils import allAmps
from database.SensorDb import SensorDb

class SensorGains(dict):
    def __init__(self, gains=None, vendorId=None, vendor='e2v',
                 db_credentials=None):
        dict.__init__(self)
        if gains is None:
            db = SensorDb(db_credentials)
            sensor = db.getSensor(vendor, vendorId)
            [self.__setitem__(amp, sensor.get_seg_result(amp, 'gain'))
             for amp in allAmps]
        else:
            try:
                [self.__setitem__(amp, gain) 
                 for amp, gain in zip(allAmps, gains)]
            except TypeError:
                [self.__setitem__(amp, gains) for amp in allAmps]

if __name__ == '__main__':
    import numpy.random as random
    # Use a single value for all segments
    g0 = SensorGains(5.5)
    for amp in allAmps:
        assert(g0[amp] == 5.5)

    # Use a list or tuple of values
    my_gains = [random.normal(5.5, 0.1) for amp in allAmps]
    g1 = SensorGains(my_gains)
    for amp in allAmps:
        assert(my_gains[amp-1] == g1[amp])

    # Check for failure.
    try:
        g2 = SensorGains()
    except KeyError, message:
        assert(("%s" % message) == "'DB_CREDENTIALS'")
        pass

    # Query the database.
    os.environ['DB_CREDENTIALS'] = '/nfs/farm/g/lsst/u1/testData/SIMData/pipeline/db_test_ro.par'
    g3 = SensorGains(vendorId='000-00', vendor='e2v')
    for amp in allAmps:
        print g3[amp]
