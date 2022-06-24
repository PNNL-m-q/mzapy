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

.. autoclass:: mzapy.calibration.LinearRegression

.. autofunction:: mzapy.calibration.LinearRegression.__init__

.. autofunction:: mzapy.calibration.LinearRegression.fit_function

.. code-block:: python3
    :caption: LinearRegression demo

    # construct some arbitrary linear data
    # slope = 3, intercept = 5, with a bit of random noise added
    import numpy as np
    X = np.arange(1, 10.1, 0.5)
    y = X * 3 + 5 + np.random.normal(size=X.shape)

    # setup linear regression
    from mzapy.calibration import LinearRegression
    lr = LinearRegression(X, y)

    # plot input values (X, y) and fitted trendline
    from matplotlib import pyplot as plt
    plt.plot(X, y, 'ko', ms=6)
    plt.plot(X, lr.transform(X), 'b-')
    plt.show()
    

