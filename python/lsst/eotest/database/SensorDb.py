"""
@brief Classes for entering sensor test results into database tables.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import lsst.eotest.image_utils as imutils
from MySQL_Database import Database

_default_callback = lambda curs : [x[0] for x in curs][0]

class NullDbObject(object):
    def __init__(self, ccdId=None, sensorDb=None):
        self.ccdId = ccdId
        self.db = sensorDb
    def add_ccd_result(self, column, value):
        pass
    def add_seg_result(self, amp, column, value):
        pass

class Sensor(object):
    def __init__(self, ccdId, sensorDb):
        self.ccdId = ccdId
        self.db = sensorDb
    def add_ccd_result(self, column, value):
        sql = ("update CCD set %s=%s where id=%i"
               % (column, value, self.ccdId))
        self.db.apply(sql)
    def add_seg_result(self, amp, column, value):
        ccdId = self.ccdId
        channelId = imutils.channelIds[amp]
        sql = """update Segment set %(column)s=%(value)s where
                 ccdId=%(ccdId)i and channelId='%(channelId)s'""" % locals()
        self.db.apply(sql)
    def get_ccd_result(self, column):
        ccdId = self.ccdId
        sql = "select %(column)s from CCD where id=%(ccdId)i" % locals()
        return self.db.apply(sql, cursorFunc=_default_callback)
    def get_seg_result(self, amp, column):
        ccdId = self.ccdId
        channelId = imutils.channelIds[amp]
        sql = """select %(column)s from Segment where ccdId=%(ccdId)i 
                 and channelId='%(channelId)s'""" % locals()
        return self.db.apply(sql, cursorFunc=_default_callback)

class SensorDbException(Exception):
    def __init__(self, *args):
        Exception.__init__(self, *args)

class SensorDb(Database):
    def __init__(self, dbdata):
        Database.__init__(self, dbdata)
    def getSensor(self, vendor, vendorId, add=False):
        try:
            return Sensor(self.getCcdId(vendor, vendorId), self)
        except IndexError:
            if add:
                return self.addSensor(vendor, vendorId)
            raise SensorDbException("Requested sensor -- vendor=%(vendor)s, vendorId=%(vendorId)s -- not in db." % locals())
    def getCcdId(self, vendor, vendorId):
        sql = """select ccdId from CCD_VendorIds where vendor='%(vendor)s'
              and vendorId='%(vendorId)s'""" % locals()
        return self.apply(sql, cursorFunc=_default_callback)
    def addSensor(self, vendor, vendorId):
        try:
            self.getCcdId(vendor, vendorId)
            raise SensorDbException("Sensor already in db.")
        except IndexError:
            pass
        sql = "insert into CCD () values ();"
        self.apply(sql)
        sql = "select id from CCD order by id desc"
        ccdId = self.apply(sql, cursorFunc=_default_callback)
        sql = """insert into CCD_VendorIds (ccdId, vendor, vendorId) values
                 (%(ccdId)i, '%(vendor)s', '%(vendorId)s')""" % locals()
        self.apply(sql)
        for amp in imutils.allAmps:
            channelId = imutils.channelIds[amp]
            sql = """insert into Segment (channelId, ccdId) values
                     ('%(channelId)s', %(ccdId)i)""" % locals()
            self.apply(sql)
        print "Adding new sensor entry for", vendor, vendorId
        return Sensor(ccdId, self)
 
if __name__ == '__main__':
    sensorDb = SensorDb('db_test_app.par')
    vendor = 'e2v'
    vendorId = '000-01'
    my_sensor = sensorDb.getSensor(vendor, vendorId, add=True)
    my_sensor.add_ccd_result('ctiSerialMean', 5e-6)
    my_sensor.add_seg_result(1, 'ctiSerial', 3e-6)

    assert(my_sensor.get_ccd_result('ctiSerialMean') == 5e-6)
    assert(my_sensor.get_seg_result(1, 'ctiSerial') == 3e-6)
    try:
        sensorDb.getSensor(vendor, '000-00', add=False)
    except SensorDbException:
        pass
