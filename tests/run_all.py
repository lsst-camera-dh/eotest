#!/usr/bin/env python
"""
@brief Discover and run all of the unit tests.  If xmlrunner.py is
available, e.g., if the jenkins subdir is in the PYTHONPATH, then
generate the xUnit compliant xml reports in the test-reports
subdirectory.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import unittest

try:
    import xmlrunner
    runner = xmlrunner.XMLTestRunner(output='test-reports', verbose=True)
except ImportError:
    runner = unittest.TextTestRunner()

loader = unittest.TestLoader()
testsuite = loader.discover('.', pattern='*tests.py')
runner.run(testsuite)
