``mzapy.calibration``
=======================================
Module for performing calibration operations.


CCS for TWIM measurements, method based on:

- Anal. Chem. 2016, 88, 7329−7336 (calibration power function)
- Nat. Protoc. 2008, 3, 1139–1152 (dt and CCS correction) 


Module Reference
---------------------------------------
.. autoclass:: mzapy.calibration.CCSCalibrationList

.. autofunction:: mzapy.calibration.CCSCalibrationList.__init__

.. autofunction:: mzapy.calibration.CCSCalibrationList.calibrated_ccs

.. autofunction:: mzapy.calibration.CCSCalibrationList.cal_curve_figure

