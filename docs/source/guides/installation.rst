Installation
==============================

From PyPI
------------------------------
``mzapy`` is installable from the PyPI via `pip <https://pip.pypa.io/en/stable/>`_:

.. code-block::

    pip install mzapy

From Source
------------------------------
The source code can also be cloned from the `GitHub repository <https://github.com/PNNL-m-q/mzapy>`_ and installed 
using `pip <https://pip.pypa.io/en/stable/>`_. This method allows for installation of development versions other than
the stable release version (``main`` branch).

.. code-block::

    # clone latest stable release (main branch)
    git clone https://github.com/PNNL-m-q/mzapy.git
    
    # OR clone the development branch
    git clone --branch dev https://github.com/PNNL-m-q/mzapy.git

    # OR clone a specific feature development branch
    git clone --branch add_multi_dim_filter https://github.com/PNNL-m-q/mzapy.git 
    
    # install
    pip install mzapy/


Dependencies
------------------------------
* ``h5py``
* ``hdf5plugin``
* ``matplotlib``
* ``numpy``
* ``pandas``
* ``scipy``
* ``matplotlib``
