"""
@brief Script to create sensor test db tables from scratch.  Existing tables
will be dropped.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
from MySQL_Database import Database

def drop_tables(db, infile='dropTables.sql'):
    for sql in open(infile):
        try:
            db.apply(sql)
        except:
            pass
        
def create_tables(db, infile='createTables.sql'):
    sql = ''.join(open(infile).readlines())
    print sql
    db.apply(sql)

if __name__ == '__main__':
#    db = Database('db_dev_u.par')
    db = Database('db_test_u.par')
    drop_tables(db)
    create_tables(db)
