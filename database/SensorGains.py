"""
@brief Provide a standard interface for accessing segment gains, either
from user input or by querying the Sensor Test database.

@author J. Chiang <jchiang@.slac.stanford.edu>
"""
import os
from database.SensorDb import SensorDb

class SensorGains(list):
    def __init__(self, gains=None, vendorId=None, vendor='e2v', nhdu=16):
        list.__init__(self)
        if gains is None:
            db = SensorDb(os.environ['DB_CREDENTIALS'])
            sensor = db.getSensor(vendor, vendorId)
            self.extend([sensor.get_seg_result(hdu, 'gain') 
                         for hdu in range(nhdu)])
        else:
            try:
                self.extend(gains)
            except TypeError:
                self.extend([gains for hdu in range(nhdu)])

if __name__ == '__main__':
    nhdu = 16

    # Use a single value for all segments
    g0 = SensorGains(5.5)
    for i in range(nhdu):
        assert(g0[i] == 5.5)

    # Use a list or tuple of values
    my_gains = (5.1, 5.2, 5.4, 5)
    g1 = SensorGains(my_gains)
    for i, gval in enumerate(g1):
        assert(my_gains[i] == gval)

    # Check for failure.
    try:
        g2 = SensorGains()
    except KeyError, message:
        assert(("%s" % message) == "'DB_CREDENTIALS'")
        pass

    # Query the database.
    os.environ['DB_CREDENTIALS'] = '/nfs/farm/g/lsst/u1/testData/SIMData/pipeline/db_test_ro.par'
    g3 = SensorGains(vendorId='000-01', vendor='e2v')
    assert(len(g3) == nhdu)
    for gval in g3:
        print gval
