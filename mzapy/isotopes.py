"""
mzapy/isotopes.py

Dylan Ross (dylan.ross@pnnl.gov)

    module for dealing with compound isotopes and masses
"""


import math
import collections


# define monoisotopic masses for selected elements
# mass corresponds to most abundant isotope
# deuterium (D) included as well since that is a common modification for MS standards
# values from: https://www.unimod.org/masses.html
_ELEMENT_MONOISO_MASS = {
    'H':   1.007825035,
    'D':   2.014101779, 
    'C':  12.,
    'N':  14.003074,
    'O':  15.99491463,
    'Na': 22.9897677,
    'P':  30.973762,
    'S':  31.9720707,
    'K':  38.9637074,
    'Se': 79.9165196,
    'He': 4.002603254,
    'Li': 	7.016003,
}


# define the masses of C, N, O, and S heavy isotopes
_HEAVY_ISOTOPE_MASS = {
    '13C': 13.00335483,
    '15N': 15.00010897,
    '18O': 17.9991603,
    '33S': 32.971458,
    '34S': 33.967867,
}


# define relative isotope abundances for C, N, O and S isotopes
# taken from Claesen, et al. JASMS 2012
_ISOTOPE_REL_ABUNDANCE = {
    'C': {
        '12C': 0.989,
        '13C': 0.011,
    },
    'N': {
        '14N': 0.996,
        '15N': 0.004,
    },
    'O': {
        '16O': 0.998,
        '18O': 0.002,
    },
    'S': {
        '32S': 0.949,
        '33S': 0.008,
        '34S': 0.043,
    },
}


def valid_element(element):
    """
    returns a ``bool`` indicating whether a specified element (``str``) is defined
    """
    # prevents users from needing to know about _ELEMENT_MONOISO_MASS
    return element in _ELEMENT_MONOISO_MASS


class MolecularFormula(collections.UserDict):
    """
    class representing a molecular formula, acts just like a dictionary mapping elements to counts but with
    some extra utilities to make them easier to add/subtract etc.

    Attributes
    ----------
    data : ``dict(str:int)``
        underlying dict mapping elements to counts
    """

    def __init__(self, *args, **kwargs):
        """
        inits a new instance of ``MolecularFormula``.
        can either be initialized empty, initialized with a dict mapping elements to counts, or kwargs with element
        names and counts as the values:

        * ``formula = MolecularFormula()`` <- empty
        * ``formula = MolecularFormula({'C': 1, 'H': 4, 'O': 1})`` <- from dict
        * ``formula = MolecularFormula(MolecularFormula(...))`` <- from MolecularFormula
        * ``formula = MolecularFormula(C=1, H=4, O=1)`` <- from kwargs
        """
        if args:
            if type(args[0]) is dict:
                super().__init__(args[0].copy())
            elif isinstance(args[0], MolecularFormula):
                # if initialized using another MolecularFormula instance, then initialize using its underlying dict
                super().__init__(args[0].data.copy())
            else:
                msg = 'MolecularFormula: __init__: initial_data must be a dict'
                raise ValueError(msg)
        else:
            # this will either be entirely empty or it will contain all of the kwargs that were provided
            super().__init__(kwargs)

    def __add__(self, other):
        """ can add two formulas together or add a dict to this formula """
        d = self.data.copy()
        for e, c in other.items():
            if e in d:
                d[e] += c
            else: 
                d[e] = c
        return MolecularFormula(d)

    def __radd__(self, other):
        """ adding formulas is commutative, re-use __add__ """
        return self.__add__(other)

    def __iadd__(self, other):
        """ implement augmented add (+=), useful for building up a formula from different pieces """
        for e, c in other.items():
            if e in self.data:
                self.data[e] += c
            else: 
                self.data[e] = c
        return self

    def __sub__(self, other):
        """ can subtract a formula or dict from this formula, but not this formula from a dict (no rsub) """
        d = self.data.copy()
        for e, c in other.items():
            if e in d:
                d[e] -= c
            else: 
                d[e] = -c
        return MolecularFormula(d)

    def __isub__(self, other):
        """ implement augmented subtract (-=) """
        for e, c in other.items():
            if e in self.data:
                self.data[e] -= c
            else: 
                self.data[e] = -c
        return self

    def __repr__(self):
        """ same as dict but with wrapping that indicates it's not the exact same """
        return 'MolecularFormula' + super().__repr__()

    def __str__(self):
        """ 
        __repr__ still outputs simular to a dict, but __str__ will return the more familiar notation for a molecular
        formula (e.g., C12H24O2N)
        """
        s = ''
        for e, c in self.data.items():
            if e[0] in '0123456789':
                e = '[' + e + ']'  # put brackets around heavy elements 
            if c != 0:  # skip 0 count elements entirely
                if c == 1:  # no number for elements with a count of 1
                    s += e
                else:  # all other cases output element, count
                    s += '{}{}'.format(e, c)
        return s

    # TODO (Dylan Ross): Add some methods for making common modifications to a molecular formula 
    #                    like methylation (H->CH3) or hydroxylation (H->OH)?

    # TODO (Dylan Ross): Add a replace method that replaces one element with another, something like:
    #                        formula.replace('H', 'D', 7)  <- replace 7 hydrogens with deuteriums
    #                    This could be useful for accounting for labeled standards as in the above example 


# define the molecular formula changes for commonly observed MS ionization states
_ADDUCT_FORMULAS = {
    '[M]+': MolecularFormula({}),
    '[M+H]+': MolecularFormula({'H': 1}),
    '[M+Na]+': MolecularFormula({'Na': 1}),
    '[M+K]+': MolecularFormula({'K': 1}),
    '[M+2K]2+': MolecularFormula({'K': 2}),
    '[M+NH4]+': MolecularFormula({'N': 1, 'H': 4}),
    '[M+H-H2O]+': MolecularFormula({'H': -1, 'O': -1}),
    '[M-H]-': MolecularFormula({'H': -1}),
    '[M+HCOO]-': MolecularFormula({'H': 1, 'C': 1, 'O': 2}),
    '[M+CH3COO]-': MolecularFormula({'H': 3, 'C': 2, 'O': 2}),
    '[M+OAc]-': MolecularFormula({'H': 3, 'C': 2, 'O': 2}),
    '[M-2H]2-': MolecularFormula({'H': -2}),
    '[M-3H]3-': MolecularFormula({'H': -3}),
    '[M+2Na-H]+': MolecularFormula({'H': -1, 'Na': 2}),
    '[M+Li]+': MolecularFormula({'Li': 1})
}
# add protonation states from +2H to +20H
for z in range(2, 21):
    _ADDUCT_FORMULAS['[M+{z}H]{z}+'.format(z=z)] = MolecularFormula({'H': z})


# define charge states for commonly observed MS ionizations
_ADDUCT_CHARGES = {
    '[M]+': 1,
    '[M+H]+': 1,
    '[M+Na]+': 1,
    '[M+K]+': 1,
    '[M+2K]2+': 2,
    '[M+NH4]+': 1,
    '[M+H-H2O]+': 1,
    '[M-H]-': 1,
    '[M+HCOO]-': 1,
    '[M+CH3COO]-': 1,
    '[M-2H]2-': 1,
    '[M-3H]3-': 1,
    '[M+2Na-H]+': 1,
}
# add protonation states from +2H to +20H
for z in range(2, 21):
    _ADDUCT_CHARGES['[M+{z}H]{z}+'.format(z=z)] = z


def valid_ms_adduct(adduct):
    """
    returns a ``bool`` indicating whether a specified MS adduct (``str``) is defined
    """
    # prevents users from needing to know about _ADDUCT_FORMULAS
    return adduct in _ADDUCT_FORMULAS


def ms_adduct_formula(neutral_formula, adduct):
    """
    modifies as molecular formula corresponding to a specified ionization state

    Parameters
    ----------
    neutral_formula : ``dict(str:int)``
        molecular formula of input neutral species as a dictionary mapping elements (str) to their counts (int)
    adduct : ``str``
        specify the type of ion

    Returns
    -------
    ion_formula : ``mzapy.isotopes.MolecularFormula``
        molecular formula of ionized species as a dictionary mapping elements (str) to their counts (int)
    """
    neutral_formula = MolecularFormula(neutral_formula)
    if adduct not in _ADDUCT_FORMULAS:
        msg = 'ms_adduct_formula: MS adduct "{}" not recognized'
        raise ValueError(msg.format(adduct))
    adduct_formula = _ADDUCT_FORMULAS[adduct]
    return neutral_formula + adduct_formula


def monoiso_mass(formula):
    """
    caculates the monoisotopic mass (assuming only most abundant isotopes) for a molecular formula

    Parameters
    ----------
    formula : ``dict(str:int)``
        molecular formula as a dictionary mapping elements (str) to their counts (int)

    Returns
    -------
    mass : ``float``
        monoisotopic mass, accurate to 6 decimal places
    """
    mass = 0.
    for element in formula:
        if element in _ELEMENT_MONOISO_MASS:
            mass += _ELEMENT_MONOISO_MASS[element] * formula[element]
        elif element in _HEAVY_ISOTOPE_MASS:
            mass += _HEAVY_ISOTOPE_MASS[element] * formula[element]
        else:
            msg = 'monoiso_mass: element "{}" not recognized'
            raise ValueError(msg.format(element))       
    # round the result just in case any funky floating point stuff happend with the multiplications
    return round(mass, 6)


def ms_adduct_mz(neutral_formula, adduct):
    """
    modifies as molecular formula corresponding to a specified ionization state, then computes m/z

    Parameters
    ----------
    neutral_formula : ``dict(str:int)``
        molecular formula of input neutral species as a dictionary mapping elements (str) to their counts (int)
    adduct : ``str``
        specify the type of ion

    Returns
    -------
    mz : ``float``
        mass to charge ratio for specified ionization state
    """
    formula = MolecularFormula(neutral_formula)
    if adduct not in _ADDUCT_FORMULAS:
        msg = 'ms_adduct_mz: MS adduct "{}" not recognized'
        raise ValueError(msg.format(adduct))
    formula += _ADDUCT_FORMULAS[adduct]
    z = _ADDUCT_CHARGES[adduct]
    return monoiso_mass(formula) / z


def _abun(element, n):
    """ returns the abundance of most abundant isotope of an element raised to coefficient n """
    return _ISOTOPE_REL_ABUNDANCE[element][{'C': '12C', 'N': '14N', 'O': '16O', 'S': '32S'}[element]]**n if n > 0 else 1.


def _abun_M1(w, x, y, z):
    """ computes the abundance for M+1 isotopes """
    a, b, c = 0., 0., 0.
    if w > 0:
        # (1x) 13C
        a = w * _ISOTOPE_REL_ABUNDANCE['C']['13C'] * _abun('C', w - 1) * _abun('N', x) * _abun('O', y) * _abun('S', z) 
    if x > 0:
        # (1x) 15N
        b = x * _ISOTOPE_REL_ABUNDANCE['N']['15N'] * _abun('C', w) * _abun('N', x - 1) * _abun('O', y) * _abun('S', z)
    if z > 0:
        # (1x) 33S
        c = z * _ISOTOPE_REL_ABUNDANCE['S']['33S'] * _abun('C', w) * _abun('N', x) * _abun('O', y) * _abun('S', z - 1)
    return a, b, c


def _mass_M1(formula, w, x, y, z):
    """ compute the mass for M+1 isotopes """
    a, b, c = 0., 0., 0.
    if w > 0:
        # (1x) 13C
        _formula = formula.copy()
        _formula['C'] -= 1
        a = _HEAVY_ISOTOPE_MASS['13C'] + monoiso_mass(_formula)
    if x > 0:
        # (1x) 15N
        _formula = formula.copy()
        _formula['N'] -= 1
        b = _HEAVY_ISOTOPE_MASS['15N'] + monoiso_mass(_formula)
    if z > 0:
        # (1x) 33S
        _formula = formula.copy()
        _formula['S'] -= 1
        c = _HEAVY_ISOTOPE_MASS['33S'] + monoiso_mass(_formula)
    return a, b, c


def _abun_M2(w, x, y, z):
    """ computes the abundance for M+2 isotopes """
    a, b, c, d, e, f, g, h = 0., 0., 0., 0., 0., 0., 0., 0.
    if w > 1:
        # (2x) 13C
        n = math.factorial(w) / (2 * math.factorial(w - 2))
        a = n * _ISOTOPE_REL_ABUNDANCE['C']['13C']**2 * _abun('C', w - 2) * _abun('N', x) * _abun('O', y) * _abun('S', z) 
    if x > 1:
        # (2x) 15N
        n = math.factorial(x) / (2 * math.factorial(x - 2))
        b = n * _ISOTOPE_REL_ABUNDANCE['N']['15N']**2 * _abun('C', w) * _abun('N', x - 2) * _abun('O', y) * _abun('S', z)
    if y > 0:
        # (1x) 18O
        c = y * _ISOTOPE_REL_ABUNDANCE['O']['18O'] * _abun('C', w) * _abun('N', x) * _abun('O', y - 1) * _abun('S', z)
    if z > 0:
        if z > 1:
            # (2x) 33S
            n = math.factorial(z) / (2 * math.factorial(z - 2))
            d = n * _ISOTOPE_REL_ABUNDANCE['S']['33S']**2 * _abun('C', w) * _abun('N', x) * _abun('O', y) * _abun('S', z - 2)
        # (1x) 34S
        e = z * _ISOTOPE_REL_ABUNDANCE['S']['34S'] * _abun('C', w) * _abun('N', x) * _abun('O', y) * _abun('S', z - 1)
    if w > 0:
        if x > 0:
            # (1x) 13C, (1x) 15N
            f = w * x * _ISOTOPE_REL_ABUNDANCE['C']['13C'] * _ISOTOPE_REL_ABUNDANCE['N']['15N'] * _abun('C', w - 1) * _abun('N', x - 1) * _abun('O', y) * _abun('S', z)
        if z > 0:
            # (1x) 13C, (1x) 33S
            g = w * z * _ISOTOPE_REL_ABUNDANCE['C']['13C'] * _ISOTOPE_REL_ABUNDANCE['S']['33S'] * _abun('C', w - 1) * _abun('N', x) * _abun('O', y) * _abun('S', z - 1)
    if x > 0 and z > 0:
        # (1x) 15N (1x) 33S
        h = x * z * _ISOTOPE_REL_ABUNDANCE['S']['33S'] * _ISOTOPE_REL_ABUNDANCE['N']['15N'] * _abun('C', w) * _abun('N', x - 1) * _abun('O', y) * _abun('S', z - 1)
    return a, b, c, d, e, f, g, h


def _mass_M2(formula, w, x, y, z):
    """ computes the mass for M+2 isotopes """
    a, b, c, d, e, f, g, h = 0., 0., 0., 0., 0., 0., 0., 0.
    if w > 1:
        # (2x) 13C
        _formula = formula.copy()
        _formula['C'] -= 2
        a = 2. * _HEAVY_ISOTOPE_MASS['13C'] + monoiso_mass(_formula)
    if x > 1:
        # (2x) 15N
        _formula = formula.copy()
        _formula['N'] -= 2
        b = 2. * _HEAVY_ISOTOPE_MASS['15N'] + monoiso_mass(_formula)
    if y > 0:
        # (1x) 18O
        _formula = formula.copy()
        _formula['O'] -= 1
        c = _HEAVY_ISOTOPE_MASS['18O'] + monoiso_mass(_formula)
    if z > 0:
        if z > 1:
            # (2x) 33S
            _formula = formula.copy()
            _formula['S'] -= 2
            d = 2. * _HEAVY_ISOTOPE_MASS['33S'] + monoiso_mass(_formula)
        # (1x) 34S
        _formula = formula.copy()
        _formula['S'] -= 1
        e = _HEAVY_ISOTOPE_MASS['34S'] + monoiso_mass(_formula)
    if w > 0:
        if x > 0:
            # (1x) 13C, (1x) 15N
            _formula = formula.copy()
            _formula['C'] -= 1
            _formula['N'] -= 1
            f = _HEAVY_ISOTOPE_MASS['13C'] + _HEAVY_ISOTOPE_MASS['15N'] + monoiso_mass(_formula)
        if z > 0:
            # (1x) 13C, (1x) 33S
            _formula = formula.copy()
            _formula['C'] -= 1
            _formula['S'] -= 1
            g = _HEAVY_ISOTOPE_MASS['13C'] + _HEAVY_ISOTOPE_MASS['33S'] + monoiso_mass(_formula)
    if x > 0 and z > 0:
        # (1x) 15N (1x) 33S
        _formula = formula.copy()
        _formula['N'] -= 1
        _formula['S'] -= 1
        h = _HEAVY_ISOTOPE_MASS['15N'] + _HEAVY_ISOTOPE_MASS['33S'] + monoiso_mass(_formula)
    return a, b, c, d, e, f, g, h


def predict_m_m1_m2(formula, relative_abundance=True):
    """
    predicts the mass and abundance (relative to M) of M, M+1, and M+2 isotopes

    isotope abundances and masses are determined using multinomial expansion but subject to 
    the following simplifying constraints:

    * only heavy isotopes 13C, 15N, 18O, 33S, and 34S are considered
    * only M, M+1, and M+2 isotope abundances are computed

    Parameters
    ----------
    formula : ``dict(str:int)``
        molecular formula as a dictionary mapping elements (str) to their counts (int)
    relative_abundance : ``bool``, default=True
        normalize the isotope abundances relative to the M isotope

    Returns
    -------
    masses : ``list(float)``
        masses of M, M+1, and M+2 isotopes
    abundances : ``list(float)``
        abundances of M, M+1, and M+2 isotopes (relative to M if relative_abundance is True)
    """
    # w, x, y, z are coefficients for C, N, O, and S from molecular formula
    w = formula['C'] if 'C' in formula else 0
    x = formula['N'] if 'N' in formula else 0
    y = formula['O'] if 'O' in formula else 0
    z = formula['S'] if 'S' in formula else 0

    # M abundance, mass
    abun_M = _abun('C', w) * _abun('N', x) * _abun('O', y) * _abun('S', z)
    mass_M = monoiso_mass(formula)

    # M+1 abundance, mass (3 possible components, a-c)
    abun_M1abc = _abun_M1(w, x, y, z)  
    abun_M1 = sum(abun_M1abc)
    mass_M1abc = _mass_M1(formula, w, x, y, z)
    mass_M1 = round(sum([_abun * _mass for _abun, _mass in zip(abun_M1abc, mass_M1abc)]) / abun_M1, 6)

    # M+2 abundance (8 possible components, a-h)
    abun_M2abcdefgh = _abun_M2(w, x, y, z)  
    abun_M2 = sum(abun_M2abcdefgh)
    mass_M2abcdefgh = _mass_M2(formula, w, x, y, z)
    mass_M2 = round(sum([_abun * _mass for _abun, _mass in zip(abun_M2abcdefgh, mass_M2abcdefgh)]) / abun_M2, 6)

    if relative_abundance:
        return (mass_M, mass_M1, mass_M2), (1., abun_M1 / abun_M, abun_M2 / abun_M)
    else:
        return (mass_M, mass_M1, mass_M2), (abun_M, abun_M1, abun_M2)

