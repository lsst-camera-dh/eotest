"""
Unit tests for bad_column code.
"""
import unittest
import numpy as np
from lsst.eotest.image_utils import bad_column

def make_column(pixel_counts):
    """
    Concatenate alternating sequences of 0's and 1's with lengths
    given by the entries in pixel_counts.  Return a list with the
    y-coordinates of the masked (value 1) pixels.
    """
    column = []
    for i, count in enumerate(pixel_counts):
        column += [i % 2]*count
    column = np.array(column)
    return list(np.where(column != 0)[0])

class BadColumnTestCase(unittest.TestCase):
    """
    TestCase class for bad_column function.
    """
    def setUp(self):
        "Create bad and ok columns of masked pixels."
        bad_pixel_counts = ((40, 30, 500),
                            (331, 20, 439),
                            (144, 60, 3421),
                            (10, 21, 309),
                            (0, 3000, 0),
                            (92, 10, 400, 33))
        self.bad_columns = (make_column(x) for x in bad_pixel_counts)
        ok_pixel_counts = ((50, 10, 40, 10, 500),
                           (0,),
                           (400, 19, 30, 19, 494),
                           (19, 19, 1, 19, 333))
        self.ok_columns = (make_column(x) for x in ok_pixel_counts)

    def tearDown(self):
        "Nothing to tear down."
        pass

    def test_bad_column(self):
        "Test the output of the bad_column function for the different cases."
        threshold = 20
        for column in self.bad_columns:
            self.assertTrue(bad_column(column, threshold))

        for column in self.ok_columns:
            self.assertFalse(bad_column(column, threshold))

if __name__ == '__main__':
    unittest.main()
