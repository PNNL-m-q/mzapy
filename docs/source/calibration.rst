``mzapy.calibration``
=======================================
Module for performing calibration operations.


CCS for single field DTIMS, method based on:

- *Anal. Chem. 2017, 89, 9048-9055*

CCS for TWIMS measurements, method based on:

- *Nat. Protoc. 2008, 3, 1139-1152* (dt and CCS correction)
- *Anal. Chem. 2016, 88, 7329-7336* (calibration power function)
- *Anal. Chem. 2020, 92, 14976-14982* (CCS calibration for SLIM)


Module Reference
---------------------------------------

Base Object
***************************************

.. autoclass :: mzapy.calibration._CalibrationBase

.. autofunction :: mzapy.calibration._CalibrationBase.fit

.. autofunction :: mzapy.calibration._CalibrationBase.transform


Mass Calibration
***************************************

.. autoclass :: mzapy.calibration.MassCalibration

.. autofunction :: mzapy.calibration.MassCalibration.__init__

.. autofunction :: mzapy.calibration.MassCalibration.calibrated_mass

.. autofunction :: mzapy.calibration.mass_calibration_from_params



DTIMS CCS Calibration (Single-Field)
***************************************

.. autoclass :: mzapy.calibration.CCSCalibrationDTsf

.. autofunction :: mzapy.calibration.CCSCalibrationDTsf.__init__

.. autofunction :: mzapy.calibration.CCSCalibrationDTsf.calibrated_ccs

.. autofunction :: mzapy.calibration.dtsf_ccs_calibration_from_params


TWIMS CCS Calibration
***************************************

.. autoclass :: mzapy.calibration.CCSCalibrationTW

.. autofunction :: mzapy.calibration.CCSCalibrationTW.__init__

.. autofunction :: mzapy.calibration.CCSCalibrationTW.calibrated_ccs

.. autofunction :: mzapy.calibration.tw_ccs_calibration_from_params
