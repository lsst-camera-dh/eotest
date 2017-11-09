"""
@brief Database interface class using MySQL.

@author J. Chiang
"""
from __future__ import absolute_import
from builtins import zip
from builtins import object
import os
import MySQLdb
from .Parfile import Parfile

db2_imp = MySQLdb

mysql_db_data = 'mysql_db_data.par'


def nullFunc(*args):
    return None


class Database(object):
    def __init__(self, dbdata=mysql_db_data):
        self.pars = Parfile(dbdata)

    def apply(self, sql, args=None, cursorFunc=nullFunc):
        sql = sql.replace('?', '%s')
        my_connection = MySQLdb.connect(**self.pars)
        cursor = my_connection.cursor()
        try:
            if args is None:
                cursor.execute(sql)
            else:
                cursor.execute(sql, args)
            results = cursorFunc(cursor)
        except db2_imp.DatabaseError as message:
            cursor.close()
            my_connection.close()
            raise db2_imp.DatabaseError(message)
        cursor.close()
        if cursorFunc is nullFunc:
            my_connection.commit()
        my_connection.close()
        return results


def getDbObjects(tablename, dbdata=mysql_db_data):
    """Return a list of entries for the specified db table"""
    sql = "SELECT * from %s" % tablename
    db = Database(dbdata)

    def cursorFunc(cursor):
        cols = [column[0] for column in cursor.description]
        entries = []
        for item in cursor:
            entries.append(dict(list(zip(cols, item))))
        return entries
    return db.apply(sql, cursorFunc=cursorFunc)


if __name__ == '__main__':
    pass
