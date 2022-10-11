``mzapy.MZA``
==============================================
This object serves as the primary interface for interacting with raw data in the MZA format.


Module Reference
---------------------------------------
.. autoclass:: mzapy.MZA

.. autofunction:: mzapy.MZA.__init__

.. note::

    The ``cache_metadata_headers`` kwarg is used to control which metadata headers are cached in memory for faster access. 
    Caching more metadata headers reduces access time, but there is a tradeoff with the memory footprint of the ``MZA`` 
    instance.  

.. code-block:: python3
    :caption: controlling which headers are cached at initialization

    from mzapy import MZA

    # load some data from ./data/example.h5, cache the default metadata headers
    h5_cache_std = MZA('./data/example.h5')

    # load some data from ./data/example.h5, do not cache any metadata headers (smallest memory footprint)
    h5_cache_none = MZA('./data/example.h5', cache_metadata_headers=[])

    # load some data from ./data/example.h5, cache all metadata headers (largest memory footprint)
    h5_cache_all = MZA('./data/example.h5', cache_metadata_headers='all')

.. autofunction:: mzapy.MZA.close

.. autofunction:: mzapy.MZA.collect_atd_arrays_by_rt_mz

.. autofunction:: mzapy.MZA.collect_dtmz_arrays_by_rt

.. autofunction:: mzapy.MZA.collect_ms1_arrays_by_rt

.. autofunction:: mzapy.MZA.collect_ms1_arrays_by_dt

.. autofunction:: mzapy.MZA.collect_ms1_arrays_by_rt_dt

.. autofunction:: mzapy.MZA.collect_ms1_df_by_rt
    
.. autofunction:: mzapy.MZA.collect_dda_ms2_df_by_rt

.. autofunction:: mzapy.MZA.collect_ms1_df_by_rt_dt

.. autofunction:: mzapy.MZA.collect_rtdt_arrays_by_mz

.. autofunction:: mzapy.MZA.collect_rtmz_arrays

.. autofunction:: mzapy.MZA.collect_rtmz_arrays_by_dt

.. autofunction:: mzapy.MZA.collect_xic_arrays_by_mz

.. autofunction:: mzapy.MZA.collect_xic_arrays_by_mz_dt

.. autofunction:: mzapy.MZA.load_scan_cache

.. autofunction:: mzapy.MZA.metadata

.. autofunction:: mzapy.MZA.save_scan_cache

