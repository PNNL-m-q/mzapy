``mzapy.isotopes``
=======================================
This module defines functions for dealing with compound isotopes and masses


The ``MolecularFormula`` class
--------------------------------------
``MolecularFormula`` is a subclass of ``collections.UserDict`` which behaves like a typical ``dict``, but supports extended functionality that is specific to molecular formulas, like addition/subtraction operations and different methods of initialization. Because ``MolecularFormula`` has the same interface as a normal ``dict``, any code which uses the old representation of molecular formulas (a plain ``dict`` mapping elements (``str``) to their counts (``int``)) can be updated to use this class instead without breaking anything. Using this class makes most of the common molecular formula operations much cleaner and simpler to implement. The ``MolecularFormula`` class also has ``__repr__`` and ``__str__`` methods implemented: 

* ``__repr__`` outputs similar to a normal ``dict`` but with "MolecularFormula" prepended to it
* ``__str__`` outputs the formula in the familiar element/count format (*e.g.*, "C4H9O2N")


A ``MolecularFormula`` may be initialized in one of four ways: 

* empty - start without any elements 
* from ``dict`` - using the previous style of molecular formula (``dict(str:int)``)
* from ``MolecularFormula`` - initialize using another ``MolecularFormula`` instance (copy the data)
* from ``kwargs`` - initialize with ``kwargs`` where the names are the elements and values are the counts

.. code-block:: python3
    :caption: Initialization

    from mzapy.isotopes import MolecularFormula

    ### empty
    formula = MolecularFormula()
    # formula.__repr__ -> 'MolecularFormula{}'
    # formula.__str__ -> ''

    ### from dict
    old_style_formula = {'C': 3, 'H': 8, 'O': 2}
    formula = MolecularFormula(old_style_formula)
    # formula.__repr__ -> 'MolecularFormula{'C': 3, 'H': 8, 'O': 2}'
    # formula.__str__ -> 'C3H8O2'

    ### from MolecularFormula
    formula = MolecularFormula({'C': 3, 'H': 8, 'O': 2})
    new_formula = MolecularFormula(formula)
    # formula.__repr__ -> 'MolecularFormula{'C': 3, 'H': 8, 'O': 2}'
    # formula.__str__ -> 'C3H8O2'

    ### from kwargs
    formula = MolecularFormula(C=3, H=8, O=2)
    # formula.__repr__ -> 'MolecularFormula{'C': 3, 'H': 8, 'O': 2}'
    # formula.__str__ -> 'C3H8O2'


``MolecularFormula`` objects support direct addition/subtraction operations with other ``MolecularFormula`` instances and also previous style of molecular formulas (``dict(str:int)``) in most cases. Addition/subtraction operations are performed element-wise and always return a ``MolecularFormula`` instance. 

.. note:: 
    
    Addition is commutative, so adding a ``MolecularFormula`` and ``dict(str:int)`` works the same in either order. This is not the case for subtraction, however, ``MolecularFormula`` - ``dict(str:int)`` works but ``dict(str:int)`` - ``MolecularFormula`` does not.


.. code-block:: python3
    :caption: Addition and Subtraction Examples

    from mzapy.isotopes import MolecularFormula

    # C4H10O2 (butyric acid)
    butyric_acid = MolecularFormula(C=4, H=8, O=2)
    # butyric_acid.__repr__ -> "MolecularFormula{'C': 4, 'H': 8, 'O': 2}"
    # butyric_acid.__str__ -> "C4H8O2"

    # deprotonate (butyrate)
    butyrate = butyric_acid - {'H': 1}
    # butyrate.__repr__ -> "MolecularFormula{'C': 4, 'H': 7, 'O': 2}"
    # butyrate.__str__ -> "C4H7O2"

    # add ammonium counter-ion
    ammonium = MolecularFormula(N=1, H=4)
    ammonium_butyrate = ammonium + butyrate
    # ammonium_butyrate.__repr__ -> "MolecularFormula{'C': 4, 'H': 11, 'O': 2, 'N'}"
    # ammonium_butyrate.__str__ -> "C4H11O2N"

    # build a hydrocarbon formula from methylene (CH2) units in a for-loop
    octane = MolecularFormula(H=1)  # start with one terminal H
    for i in range(8):
        octane += {'C': 1, 'H': 2}  # add methylene units
    octane += {'H': 1}  # finish off with the other terminal H
    # octane.__repr__ -> "MolecularFormula{'C': 8, 'H': 18}"
    # octane.__str__ -> "C8H18"


Elements
--------
The following table summarizes the currently defined elements in the ``mzapy.isotopes`` module, 
along with the exact masses of their most abundant isotope (source: https://www.unimod.org/masses.html)

============  ================
 Element      Exact Mass   
============  ================
H               1.007825035
D               2.014101779 
C               12.0000000
N               14.003074
O               15.99491463
Na              22.9897677
P               30.973762
S               31.9720707
K               38.9637074
Se              79.9165196
He              4.002603254
Li 	            7.016003
============  ================


Additional elements may be defined by adding entries to the ``mzapy.isotopes._ELEMENT_MONOISO_MASS`` 
dictionary.


MS Adducts
----------
Molecular formulas and m/z values can be computed for various MS adducts using :func:`mzapy.isotopes.ms_adduct_formula` and 
:func:`mzapy.isotopes.ms_adduct_mz`. The available adduct types are:


============  ======== 
 adduct         z   
============  ======== 
[M]+            1
[M+H]+          1
[M+Na]+         1
[M+K]+          1
[M+2K]2+        2
[M+NH4]+        1
[M+H-H2O]+      1
[M-H]-          1
[M+HCOO]-       1
[M+CH3COO]-     1
[M-2H]2-        1
[M-3H]3-        1
[M+2Na-H]+      1
[M+2H]2+        2
[M+3H]3+        3
[M+4H]4+        4
[M+5H]5+        5
[M+6H]6+        6
[M+7H]7+        7
[M+8H]8+        8
[M+9H]9+        9
[M+10H]10+      10
[M+11H]11+      11
[M+12H]12+      12
[M+13H]13+      13
[M+14H]14+      14
[M+15H]15+      15
[M+16H]16+      16
[M+17H]17+      17
[M+18H]18+      18
[M+19H]19+      19
[M+20H]20+      20
============  ======== 


Module Reference
---------------------------------------

Molecular Formula Object
***************************************

.. autoclass :: mzapy.isotopes.MolecularFormula

.. autofunction :: mzapy.isotopes.MolecularFormula.__init__

Utility Functions
***************************************

.. autofunction :: mzapy.isotopes.valid_element

.. autofunction :: mzapy.isotopes.valid_ms_adduct

.. autofunction :: mzapy.isotopes.monoiso_mass

.. autofunction :: mzapy.isotopes.ms_adduct_formula

.. autofunction :: mzapy.isotopes.ms_adduct_mz

.. autofunction :: mzapy.isotopes.predict_m_m1_m2

