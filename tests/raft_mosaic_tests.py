"""
Test code for raft-level mosaicking code.
"""
from __future__ import absolute_import, print_function
import os
import glob
import unittest
import lsst.eotest.raft as raftTest

class RaftMosaicTestCase(unittest.TestCase):
    "Test case class for raft_image_mosaic module."
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_raft_image_mosaic(self):
        root_dir = '/nfs/farm/g/lsst/u1/jobHarness/jh_archive-test/LCA-10753_RSA/LCA-10753_RSA-003_ETU1-Dev/4578D/qe_raft_acq_sim/v0/25597'
        infiles = glob.glob(os.path.join(root_dir, 'S??',
                                         '*_lambda_flat_1000_*.fits'))
        raft_mosaic = raftTest.RaftMosaic(infiles)

if __name__ == '__main__':
    unittest.main()
