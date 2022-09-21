"""
mzapy/test/peaks.py

Dylan Ross (dylan.ross@pnnl.gov)

    unit tests for peaks module
"""


import unittest

import numpy as np

from mzapy.peaks import lerp_1d, lerp_2d


class TestSignalProcessing(unittest.TestCase):
    """ unit tests for signal processing functions (interpolation, smoothing, etc.) """

    def test_lerp_1d_bad_xy_len(self):
        # test cases where x/y data are the wrong lengths
        with self.assertRaises(ValueError):
            # empty x and y arrays
            x, y = np.array([]), np.array([])
            lerp_1d(x, y, 0., 1., 100)
        with self.assertRaises(ValueError):
            # x and y arrays too short
            x, y = np.array([1]), np.array([1])
            lerp_1d(x, y, 0., 1., 100)
        with self.assertRaises(ValueError):
            # x and y arrays have different lengths
            x, y = np.array([1, 2, 3]), np.array([1, 2])
            lerp_1d(x, y, 0., 1., 100)


class TestPeakFitting(unittest.TestCase):
    """ unit tests for peak fitting functions """

    def test_find_peaks_1d_localmax_bad_input(self):
        # tests localmax peak fitting function with bad inputs 
        pass


if __name__ == '__main__':
    # run all TestCases
    unittest.main(verbosity=2)


