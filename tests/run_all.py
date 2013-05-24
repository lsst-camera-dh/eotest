#!/usr/bin/env python
"""
@brief Discover and run all of the unit tests.  If xmlrunner.py is
available, e.g., if the jenkins subdir is in the PYTHONPATH, then
generate the xUnit compliant xml reports in the test-reports
subdirectory.

@author J. Chiang <jchiang@slac.stanford.edu>
"""
import unittest

loader = unittest.TestLoader()
testsuite = loader.discover('.', pattern='*tests.py')
print "Found %i tests" % testsuite.countTestCases()

try:
    import xmlrunner
    runner = xmlrunner.XMLTestRunner(output='test-reports', verbose=True)
    runner.run(testsuite)
except ImportError:
    runner = unittest.TextTestRunner()
    runner.run(testsuite)
