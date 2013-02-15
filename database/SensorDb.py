"""
@brief Classes for entering sensor test results into database tables.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
from MySQL_Database import Database

class NullDbObject(object):
    def __init__(self, ccdId, sensorDb):
        self.ccdId = ccdId
        self.db = sensorDb
    def add_ccd_result(self, column, value):
        pass
    def add_seg_result(self, segment, column, value):
        pass

class Sensor(object):
    def __init__(self, ccdId, sensorDb):
        self.ccdId = ccdId
        self.db = sensorDb
    def add_ccd_result(self, column, value):
        sql = ("update CCD set %s=%s where id=%i"
               % (column, value, self.ccdId))
        self.db.apply(sql)
    def add_seg_result(self, segment, column, value):
        sql = ("""update Segment set %s=%s where
               ccdId=%i and channelId='%02o'"""
               % (column, value, self.ccdId, segment))
        self.db.apply(sql)

_default_callback = lambda curs : [x[0] for x in curs][0]

class SensorDbException(Exception):
    def __init__(self, *args):
        Exception.__init__(self, *args)

class SensorDb(Database):
    def __init__(self, dbdata):
        Database.__init__(self, dbdata)
    def getSensor(self, vendor, vendorId):
        try:
            return self._addSensor(vendor, vendorId)
        except SensorDbException:
            sql = """select ccdId from CCD_VendorIds
                  where vendor='%(vendor)s' and
                  vendorId='%(vendorId)s'""" % locals()
            return Sensor(self.apply(sql, cursorFunc=_default_callback), self)
    def _addSensor(self, vendor, vendorId):
        sql = """select ccdId from CCD_VendorIds where vendor='%(vendor)s'
              and vendorId='%(vendorId)s'""" % locals()
        try:
            ccdId = self.apply(sql, cursorFunc=_default_callback)
            raise SensorDbException("Sensor already in db")
        except IndexError:
            pass
        sql = "insert into CCD () values ();"
        self.apply(sql)
        sql = "select id from CCD order by id desc"
        ccdId = self.apply(sql, cursorFunc=_default_callback)
        sql = """insert into CCD_VendorIds (ccdId, vendor, vendorId) values
                 (%(ccdId)i, '%(vendor)s', '%(vendorId)s')""" % locals()
        self.apply(sql)
        for segment in range(16):
            channelId = "%02o" % segment
            sql = """insert into Segment (channelId, ccdId) values
                     ('%(channelId)s', %(ccdId)i)""" % locals()
            self.apply(sql)
        print "Adding new sensor entry for", vendor, vendorId
        return Sensor(ccdId, self)
 
if __name__ == '__main__':
    sensorDb = SensorDb('mysql_db_data_app.par')
    vendor = 'e2v'
    vendorId = '000-01'
    my_sensor = sensorDb.getSensor(vendor, vendorId)
    my_sensor.add_ccd_result('ctiSerialMean', 5e-6)
    my_sensor.add_seg_result(0, 'ctiSerial', 3e-6)
