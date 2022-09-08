"""
mzapy/test/isotope.py

Dylan Ross (dylan.ross@pnnl.gov)

    unit tests for isotopes module
"""


import unittest

from mzapy.isotopes import MolecularFormula


class TestMolecularFormula(unittest.TestCase):
    """ 
    unit tests for the MolecularFormula class 
    
    test groups:
        0 - initialization
        1 - add/subtract operations
        2 - repr/str methods
        3 - modifying MolecularFormula does not affect other MolecularFormula
    """

    def test_0_init_empty(self):
        mol = MolecularFormula()
        self.assertDictEqual(mol.data, {})

    def test_0_init_from_dict(self):
        d = {'C': 1, 'H': 4, 'O': 1}
        mol = MolecularFormula(d)
        self.assertDictEqual(mol.data, {'C': 1, 'H': 4, 'O': 1})

    def test_0_init_from_molecular_formula(self):
        mol1 = MolecularFormula({'C': 1, 'H': 4, 'O': 1})
        mol2 = MolecularFormula(mol1)
        self.assertDictEqual(mol1.data, mol2.data)

    def test_0_init_from_kwargs(self):
        mol = MolecularFormula(C=1, H=4, O=1)
        self.assertDictEqual(mol.data, {'C': 1, 'H': 4, 'O': 1})

    def test_1_add_molecular_formulas(self):
        mol1 = MolecularFormula(C=1, H=4, O=1)
        mol2 = MolecularFormula(C=2, H=4)
        mol3 = mol1 + mol2
        self.assertDictEqual(mol3.data, {'C': 3, 'H': 8, 'O': 1})

    def test_1_add_dict_to_molecular_formula(self):
        mol = MolecularFormula(C=1, H=4, O=1)
        d = {'C': 2, 'H': 4}
        mol3 = mol + d
        self.assertDictEqual(mol3.data, {'C': 3, 'H': 8, 'O': 1})

    def test_1_radd_dict_to_molecular_formula(self):
        mol = MolecularFormula(C=1, H=4, O=1)
        d = {'C': 2, 'H': 4}
        mol3 = d + mol
        self.assertDictEqual(mol3.data, {'C': 3, 'H': 8, 'O': 1})

    def test_1_iadd_molecular_formulas(self):
        mol1 = MolecularFormula(C=1, H=4, O=1)
        mol2 = MolecularFormula(C=2, H=4)
        mol1 += mol2
        self.assertDictEqual(mol1.data, {'C': 3, 'H': 8, 'O': 1})

    def test_1_iadd_dict_to_molecular_formula(self):
        mol = MolecularFormula(C=1, H=4, O=1)
        d = {'C': 2, 'H': 4}
        mol += d
        self.assertDictEqual(mol.data, {'C': 3, 'H': 8, 'O': 1})

    def test_1_sub_molecular_formulas(self):
        mol1 = MolecularFormula(C=3, H=8, O=1)
        mol2 = MolecularFormula(C=2, H=4)
        mol3 = mol1 - mol2
        self.assertDictEqual(mol3.data, {'C': 1, 'H': 4, 'O': 1})

    def test_1_sub_dict_from_molecular_formula(self):
        mol = MolecularFormula(C=3, H=8, O=1)
        d = {'C': 2, 'H': 4}
        mol3 = mol - d
        self.assertDictEqual(mol3.data, {'C': 1, 'H': 4, 'O': 1})
    
    def test_1_sub_molecular_formula_from_dict(self):
        mol = MolecularFormula(C=3, H=8, O=1)
        d = {'C': 2, 'H': 4}
        with self.assertRaises(TypeError):
            # not allowed to subtract MolecularFormula from dict (no __rsub__)
            mol3 = d - mol
    
    def test_1_isub_molecular_formulas(self):
        mol1 = MolecularFormula(C=3, H=8, O=1)
        mol2 = MolecularFormula(C=2, H=4)
        mol1 -= mol2
        self.assertDictEqual(mol1.data, {'C': 1, 'H': 4, 'O': 1})

    def test_1_sub_dict_from_molecular_formula(self):
        mol = MolecularFormula(C=3, H=8, O=1)
        d = {'C': 2, 'H': 4}
        mol -= d
        self.assertDictEqual(mol.data, {'C': 1, 'H': 4, 'O': 1})

    def test_2_repr(self):
        # TODO(Dylan Ross): add some more cases here
        self.assertEqual(MolecularFormula(C=4, H=12, O=1).__repr__(), "MolecularFormula{'C': 4, 'H': 12, 'O': 1}")

    def test_2_str(self):
        # TODO(Dylan Ross): add some more cases here
        self.assertEqual(MolecularFormula(C=4, H=12, O=1).__str__(), "C4H12O")
        self.assertEqual(MolecularFormula({'C': 4, 'H': 12, '18O': 1}).__str__(), "C4H12[18O]")
    
    def test_3_molecular_formulas_are_independent(self):
        d = {'C': 2, 'H': 4}
        mol = MolecularFormula(d)
        O = MolecularFormula(O=1)
        mol += O
        # mol should have had O added to it and d and O (MolecularFormula) should both be unchanged after the addition
        self.assertDictEqual(d, {'C': 2, 'H': 4})
        self.assertDictEqual(O.data, {'O': 1})
        self.assertDictEqual(mol.data, {'C': 2, 'H': 4, 'O': 1})


if __name__ == '__main__':
    # run all TestCases
    unittest.main(verbosity=2)

