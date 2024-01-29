"""
mzapy/test/calibration.py

Dylan Ross (dylan.ross@pnnl.gov)

    unit tests for calibration module
"""


import unittest

from mzapy.calibration import (
    mass_calibration_from_params,
    dtsf_ccs_calibration_from_params,
    tw_ccs_calibration_from_params,
)


class TestCalFactoryFuncs(unittest.TestCase):
    """ unit tests for the calibration class factory functions """
    
    def test_mass_cal_from_params(self):
        # should be no errors initializing with class factory
        cal = mass_calibration_from_params("linear", (1., 0.))
        # make sure the instance works as expected
        # AlmostEqual in case of any floating point stuff causing slight difference
        self.assertAlmostEqual(cal.calibrated_mass(420.), 420.,
                               msg="got wrong value from cal.calibrated_mass")

    def test_mass_cal_from_params_bad_fit_func(self):
        # should raise a ValueError 
        with self.assertRaises(ValueError,
                               msg="bad fit function should raise ValueError"):
            cal = mass_calibration_from_params("not good", (1., 0.))

    def test_dtsf_ccs_cal_from_params(self):
        # should be no errors initializing with class factory
        cal = dtsf_ccs_calibration_from_params(1, 0., 1.0)
        # need a real value to test..?
        #self.assertAlmostEqual(cal.calibrated_ccs(,), )

    def test_tw_ccs_cal_from_params(self):
        # should be no errors initializing with class factory
        cal = tw_ccs_calibration_from_params(1, "linear", (0.08, 400.))
        # make sure the instance works as expected
        # AlmostEqual in case of any floating point stuff causing slight difference
        self.assertAlmostEqual(cal.calibrated_ccs(790.4457, 15584.), 316.6311342,
                               msg="got wrong value from cal.calibrated_ccs")

    def test_tw_ccs_cal_from_params_bad_fit_func(self):
        # should raise a ValueError 
        with self.assertRaises(ValueError,
                               msg="bad fit function should raise ValueError"):
            cal = tw_ccs_calibration_from_params(1, "not good", (1., 0.))


if __name__ == '__main__':
    # run all TestCases
    unittest.main(verbosity=2)

