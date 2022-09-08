"""
mzapy/test/__main__.py

Dylan Ross (dylan.ross@pnnl.gov)

    runs unit tests
"""


import unittest

from mzapy.test.isotopes import TestMolecularFormula


if __name__ == '__main__':
    # run all imported TestCases
    unittest.main(verbosity=2)