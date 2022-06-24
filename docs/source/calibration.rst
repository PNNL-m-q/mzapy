``mzapy.calibration``
=======================================
Module for performing calibration operations.


CCS for TWIM measurements, method based on:

- *Nat. Protoc. 2008, 3, 1139-1152* (dt and CCS correction)
- *Anal. Chem. 2016, 88, 7329-7336* (calibration power function)
- *Anal. Chem. 2020, 92, 14976-14982* (CCS calibration for SLIM)


Module Reference
---------------------------------------

.. autoclass:: mzapy.calibration._CalibrationBase

.. autofunction:: mzapy.calibration._CalibrationBase.fit

.. autofunction:: mzapy.calibration._CalibrationBase.transform

.. autoclass:: mzapy.calibration.TWCCSCalibration

.. autofunction:: mzapy.calibration.TWCCSCalibration.__init__

.. autofunction:: mzapy.calibration.TWCCSCalibration.calibrated_ccs
