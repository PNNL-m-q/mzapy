``mzapy.dda``
==============================================
This module contains a couple of specialized classes for extracting data from DDA MS/MS MZA files.


Module Reference
---------------------------------------

.. note::
    The ``MsmsReaderDdaCachedMs1`` class has the same interface as the base ``MsmsReaderDda`` class
    except that it has additional attributes (``arrays_mz_cached` and ``arrays_i_cached``) which store
    cached MS1 data to make data extraction faster. Use of this class should be avioded in contexts
    with multiprocessing since the cached data can have a significant memory footprint. 


Initialization
***************************************

.. autoclass :: mzapy.dda.MsmsReaderDda

.. autoclass :: mzapy.dda.MsmsReaderDdaCachedMs1

.. autofunction :: mzapy.dda.MsmsReaderDda.__init__

.. note:: 
    The ``drop_scans`` argument takes a list of scan numbers to drop from the data file during reading. 
    This is used to get rid of erroneous scans that did not get properly recorded in the original data
    file and that interfere with extracting data. The argument is optional, and if it is ``None``
    then no scans are dropped. 

.. autofunction :: mzapy.dda.MsmsReaderDda.close


Data Extraction
***************************************

.. autofunction :: mzapy.dda.MsmsReaderDda.get_pre_mzs

.. autofunction :: mzapy.dda.MsmsReaderDda.get_tic

.. autofunction :: mzapy.dda.MsmsReaderDda.get_chrom

.. autofunction :: mzapy.dda.MsmsReaderDda.get_msms_spectrum
