"""
@brief Unit tests for Parfile.py module.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import unittest
from lsst.eotest.database.Parfile import Parfile


class ParfileTestCase(unittest.TestCase):
    def setUp(self):
        try:
            os.remove('foo.pars')
        except OSError:
            pass
        self.test_file = 'pars_test.txt'
        self.write_test_file = 'write_test.txt'
        output = open(self.test_file, 'w')
        output.write('ft1file = ft1.fits\n')
        output.close()

    def tearDown(self):
        os.remove(self.test_file)
        try:
            os.remove(self.write_test_file)
        except OSError:
            pass

    def testEmptyFile(self):
        self.assertRaises(IOError, Parfile, 'foo.pars')

    def testFixedKeys(self):
        pars = Parfile('foo.pars', fixed_keys=False)
        pars['ft2file'] = 'ft2.fits'
        pars.write()
        pars = Parfile('foo.pars')
        self.assertRaises(KeyError, pars.__setitem__, 'ft1file', '')

    def testReadWrite(self):
        pars = Parfile('foo.pars', fixed_keys=False)
        pars['ft1file'] = 'ft1.fits'
        pars['ft2file'] = 'ft2.fits'
        pars['ra'] = 83.57
        pars['dec'] = 22.01
        pars['ft1file'] = 'ft1_file.fits'
        pars.write()
        foo = Parfile('foo.pars')
        self.assertEquals(foo['ft1file'], 'ft1_file.fits')
        self.assertEquals(foo['ft2file'], 'ft2.fits')
        self.assertEquals(foo['ra'], 83.57)
        self.assertEquals(foo['dec'], 22.01)

    def testWriteMethod(self):
        outfile = self.write_test_file
        infile = self.test_file
        pars = Parfile(infile)
        pars.write(outfile)
        bar = Parfile(infile)
        foo = Parfile(outfile)
        self.assertEquals(foo, bar)
        foo['ft1file'] = 'ft2.fits'
        self.assertNotEquals(foo, bar)


if __name__ == '__main__':
    unittest.main()
