"""
@brief Unit tests for AmplifierGeometry.py module.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import os
import unittest

try:
    from lsst.eotest.sensor import AmplifierGeometry
except ImportError:
    # This is to allow this unit test to run on the inadequately
    # configured lsst-build01 on which Jenkins at SLAC runs.
    print "Error importing lsst.eotest.sensor"
    import sys
    sys.path.insert(0, os.path.join(os.environ['TEST_SCRIPTS_DIR'],
                                    'python', 'lsst', 'eotest', 'sensor'))
    from AmplifierGeometry import AmplifierGeometry

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
        self.e2v = AmplifierGeometry(prescan=10,
                                     serial_overscan=20, parallel_overscan=20,
                                     detxsize=4336, detysize=4044)
    def tearDown(self):
        pass
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
        
if __name__ == '__main__':
    unittest.main()
