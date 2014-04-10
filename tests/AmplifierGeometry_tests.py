"""
@brief Unit tests for AmplifierGeometry.py module.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import unittest

from lsst.eotest.sensor import AmplifierGeometry, makeAmplifierGeometry, amp_loc
import lsst.eotest.sensor.sim_tools as sim_tools
import lsst.afw.geom as afwGeom

e2v_amp_geom = dict(
    [(1, dict([('DATASEC', '[11:522,1:2002]'),
               ('DETSEC', '[512:1,1:2002]'),
               ('DETSIZE', '[1:4096,1:4004]')])),
     (2, dict([('DATASEC', '[11:522,1:2002]'),
               ('DETSEC', '[1024:513,1:2002]'),
               ('DETSIZE', '[1:4096,1:4004]')])),
     (3, dict([('DATASEC', '[11:522,1:2002]'),
               ('DETSEC', '[1536:1025,1:2002]'),
               ('DETSIZE', '[1:4096,1:4004]')])),
     (4, dict([('DATASEC', '[11:522,1:2002]'),
               ('DETSEC', '[2048:1537,1:2002]'),
               ('DETSIZE', '[1:4096,1:4004]')])),
     (5, dict([('DATASEC', '[11:522,1:2002]'),
               ('DETSEC', '[2560:2049,1:2002]'),
               ('DETSIZE', '[1:4096,1:4004]')])),
     (6, dict([('DATASEC', '[11:522,1:2002]'),
               ('DETSEC', '[3072:2561,1:2002]'),
               ('DETSIZE', '[1:4096,1:4004]')])),
     (7, dict([('DATASEC', '[11:522,1:2002]'),
               ('DETSEC', '[3584:3073,1:2002]'),
               ('DETSIZE', '[1:4096,1:4004]')])),
     (8, dict([('DATASEC', '[11:522,1:2002]'),
               ('DETSEC', '[4096:3585,1:2002]'),
               ('DETSIZE', '[1:4096,1:4004]')])),
     (9, dict([('DATASEC', '[11:522,1:2002]'),
               ('DETSEC', '[3585:4096,4004:2003]'),
               ('DETSIZE', '[1:4096,1:4004]')])),
     (10, dict([('DATASEC', '[11:522,1:2002]'),
                ('DETSEC', '[3073:3584,4004:2003]'),
                ('DETSIZE', '[1:4096,1:4004]')])),
     (11, dict([('DATASEC', '[11:522,1:2002]'),
                ('DETSEC', '[2561:3072,4004:2003]'),
                ('DETSIZE', '[1:4096,1:4004]')])),
     (12, dict([('DATASEC', '[11:522,1:2002]'),
                ('DETSEC', '[2049:2560,4004:2003]'),
                ('DETSIZE', '[1:4096,1:4004]')])),
     (13, dict([('DATASEC', '[11:522,1:2002]'),
                ('DETSEC', '[1537:2048,4004:2003]'),
                ('DETSIZE', '[1:4096,1:4004]')])),
     (14, dict([('DATASEC', '[11:522,1:2002]'),
                ('DETSEC', '[1025:1536,4004:2003]'),
                ('DETSIZE', '[1:4096,1:4004]')])),
     (15, dict([('DATASEC', '[11:522,1:2002]'),
                ('DETSEC', '[513:1024,4004:2003]'),
                ('DETSIZE', '[1:4096,1:4004]')])),
     (16, dict([('DATASEC', '[11:522,1:2002]'),
                ('DETSEC', '[1:512,4004:2003]'),
                ('DETSIZE', '[1:4096,1:4004]')]))]
    )

full_segment = afwGeom.Box2I(afwGeom.Point2I(0, 0),
                             afwGeom.Point2I(541, 2021))

prescan = afwGeom.Box2I(afwGeom.Point2I(0, 0),
                        afwGeom.Point2I(9, 2021))

imaging = afwGeom.Box2I(afwGeom.Point2I(10, 0),
                        afwGeom.Point2I(521, 2001))

serial_overscan = afwGeom.Box2I(afwGeom.Point2I(522, 0), 
                                afwGeom.Point2I(541, 2021))

parallel_overscan = afwGeom.Box2I(afwGeom.Point2I(10, 2002),
                                  afwGeom.Point2I(522, 2021))

class AmplifierGeometryTestCase(unittest.TestCase):
    def setUp(self):
        self.e2v = AmplifierGeometry(prescan=10, nx=512, ny=2002,
                                     detxsize=4336, detysize=4044,
                                     amp_loc=amp_loc['E2V'])
        self.e2v_test_file = 'test_e2v_image.fits'
        ccd0 = sim_tools.CCD(geometry=self.e2v)
        ccd0.writeto(self.e2v_test_file, bitpix=16)

        self.itl = AmplifierGeometry(prescan=3, nx=509, ny=2000,
                                     amp_loc=amp_loc['ITL'])
        self.itl_test_file = 'test_itl_image.fits'
        ccd1 = sim_tools.CCD(geometry=self.itl)
        ccd1.writeto(self.itl_test_file, bitpix=16)
    def tearDown(self):
        os.remove(self.itl_test_file)
        os.remove(self.e2v_test_file)
    def test_e2v_keywords(self):
        self.assertEquals(self.e2v.DETSIZE, '[1:4336,1:4044]')
        for amp in range(1, 17):
            for key in ('DATASEC', 'DETSEC', 'DETSIZE'):
                self.assertEquals(self.e2v[amp][key], e2v_amp_geom[amp][key])
    def test_e2v_geometry(self):
        self.assertEquals(self.e2v.full_segment, full_segment)
        self.assertEquals(self.e2v.prescan, prescan)
        self.assertEquals(self.e2v.imaging, imaging)
        self.assertEquals(self.e2v.serial_overscan, serial_overscan)
        self.assertEquals(self.e2v.parallel_overscan, parallel_overscan)
    def test_makeAmplifierGeometry_factory(self):
        e2v_geom = makeAmplifierGeometry(self.e2v_test_file)
        self.assertEquals(self.e2v, e2v_geom)

        itl_geom = makeAmplifierGeometry(self.itl_test_file)
        self.assertEquals(self.itl, itl_geom)

        self.assertNotEqual(self.e2v, itl_geom)
        self.assertNotEqual(self.itl, e2v_geom)
        
if __name__ == '__main__':
    unittest.main()
