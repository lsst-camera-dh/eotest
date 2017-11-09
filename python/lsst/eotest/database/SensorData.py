"""
@brief Provide access to db quantities for a specified CCD.
Properties may be plotted versus segment id.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import numpy as np
import pylab
from image_utils import allAmps, channelIds
from database.SensorDb import SensorDb


class SensorData(object):
    def __init__(self, sensor_id, vendor='e2v', dbdata=None):
        self.db = SensorDb(dbdata)
        self.sensor = self.db.getSensor(vendor, sensor_id)
        self._get_properties()
        self._get_data()

    def _get_properties(self):
        sql = """select column_name from information_schema.COLUMNS
                 where table_schema='%s'
                 and table_name='Segment'""" % self.db.pars['db']

        def func(curs): return [x[0] for x in curs]
        self.properties = self.db.apply(sql, cursorFunc=func)

    def _get_data(self):
        sql = """select * from Segment where
                 ccdId=%s order by id asc""" % self.sensor.ccdId

        def func(curs): return [x for x in curs]
        coldata = np.array(self.db.apply(sql, cursorFunc=func)).transpose()
        self.data = dict([(name, column) for name, column
                          in zip(self.properties, coldata)])
        self.indx = self.data['channelId'].argsort()
        self.nchan = len(self['channelId'])

    def __getitem__(self, property):
        "Return the column values ordered by channelId."
        return self.data[property][self.indx]

    def addColumn(self, property, values):
        self.data[property] = np.array(values)

    def plot(self, property, yrange=None, ylabel=None, title=None,
             outfile=None):
        fig = pylab.figure()
        if title is None:
            title = '%s %s' % (vendor, sensor_id)
        fig.suptitle(title, y=0.95)
        pylab.plot(range(self.nchan), self[property], 'k.')
        pylab.xticks(range(self.nchan), self['channelId'])
        pylab.xlabel('Segment')
        if ylabel is not None:
            pylab.ylabel(ylabel)
        else:
            pylab.ylabel(property)
        axisrange = list(pylab.axis())
        axisrange[:2] = (-1, self.nchan)
        if yrange is not None:
            axisrange[2:] = yrange[:2]
        pylab.axis(axisrange)
        pylab.show()
        if outfile is not None:
            pylab.savefig(outfile)
        return fig


if __name__ == '__main__':
    dbdata = 'db_test_ro.par'   # This file contains the mysql credentials.
    vendor = 'e2v'
    sensor_id = '000-00'

    sensor_data = SensorData(sensor_id, dbdata=dbdata)

    #
    # Make some plots
    #
    sensor_data.plot('gain', yrange=(0, 6), ylabel='System Gain (e-/DN)',
                     outfile='000_00_gain.png')
    sensor_data.plot('readNoise', ylabel='Read Noise (e- rms)',
                     outfile='000_00_read_noise.png')
    sensor_data.plot('maxDeviation',
                     ylabel='Maximum fractional deviation from linear',
                     outfile='000_00_max_deviation.png')
    sensor_data.plot('fullWell', ylabel='Full Well (e-)',
                     outfile='000_00_full_well.png')

    #
    # Compute a new column and add to the SensorData object
    #
    maxdev_percent = sensor_data.data['maxDeviation']*100.
    sensor_data.addColumn('maxdev_percent', maxdev_percent)
    sensor_data.plot('maxdev_percent', yrange=(0, 2),
                     ylabel='Maximum fractional deviation from linear (%)',
                     outfile='000_00_maxdev_percent.png')
