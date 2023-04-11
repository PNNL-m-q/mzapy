"""
mzapy/calibration.py
    
Dylan Ross (dylan.ross@pnnl.gov)

    Module for performing calibration operations.

    CCS for TWIM measurements, method based on:
        - Nat. Protoc. 2008, 3, 1139-1152 (dt and CCS correction)
        - Anal. Chem. 2016, 88, 7329-7336 (calibration power function)
        - Anal. Chem. 2020, 92, 14976-14982 (CCS calibration for SLIM)
"""


from scipy import optimize
import numpy as np

from mzapy.isotopes import monoiso_mass


class _CalibrationBase:
    """
    Modular base object for performing calibration

    Works similar to BaseEstimator from sklearn, implements fit and transform methods. 
    Subclasses must implement self.fit_function with signature: ``self.fit_function(X, *params) -> y``
    and they must override the attribute self.init_params with initial parameter values for the fitting function.

    Methods
    -------
    fit(X, y) -> y_fit
        Takes input values (X) and known output values (y), then performs least-squares optimization 
        with subclass defined self.fit_function. Returns the fitted values (y_fit) 
        and sets self.opt_params with the optimized parameters
    transform(X) -> y_transform
        Takes  input values (X) and uses subclass defined self.fit_function with optimized parameters 
        (self.opt_params) to produce output values (y_transform). Requires self.fit method to run 
        successfully first so that self.opt_params gets set

    Attributes
    ----------
    init_params : ``tuple(?)``
        initial parameter values for fitting function, must be set by subclass
    opt_params : ``tuple(?)``
        optimized parameter values for fitting function, set by self.fit method, initially set to None to indicate
        the calibration has not been fit yet
    """

    def __init__(self):
        """
        initialize new _CalibrationBase instance
        """
        # set attributes
        self.init_params = None
        self.opt_params = None  # should be overridden in subclasses

    @staticmethod
    def fit_function(X, *params):
        """
        ! must be overridden by subclass, calling this function raises an exception !

        Parameters
        ----------
        X : ``numpy.ndarray(float)`` or ``float``
            input values
        params : ``tuple(?)``
            fitting function parameters

        Returns
        -------
        y : ``numpy.ndarray(float)`` or ``float``
            fitting function output
        """
        msg = ('_CalibrationBase: fit_function: This is the base class fitting function,'
               ' it should have been overridden by subclass')
        raise RuntimeError(msg)

    def fit(self, X, y):
        """
        Takes input values (X) and known output values (y), then performs least-squares optimization 
        with subclass defined self.fit_function. Returns the fitted values (y_fit) 
        and sets self.opt_params with the optimized parameters

        Parameters
        ----------
        X : ``numpy.ndarray(float)``
            input values
        y : ``numpy.ndarray(float)``
            known output values

        Returns
        -------
        y_fit : ``numpy.ndarray(float)``
            fitted output values
        """
        # check that self.init_params has been set
        if self.init_params is None:
            msg = '_CalibrationBase: fit: self.init_params has not been set, this should have been set by subclass'
            raise RuntimeError(msg)
        try:
            opt_params, cov = optimize.curve_fit(self.fit_function, X, y, maxfev=50000, p0=self.init_params)
        except Exception as e:
            msg = '_CalibrationBase: fit: unable to fit calibration:\n{}'
            raise RuntimeError(msg.format(e))            
        # if fit successful, set the optimized parameters
        self.opt_params = opt_params
        # then return the transformed input values (y_fit)
        return self.transform(X)

    def transform(self, X):
        """
        Takes  input values (X) and uses subclass defined self.fit_function with optimized parameters 
        (self.opt_params) to produce output values (y_transform). Requires self.fit method to run 
        successfully first so that self.opt_params gets set

        Parameters
        ----------
        X : ``numpy.ndarray(float)`` or ``float``
            input values

        Returns
        -------
        y_transform : ``numpy.ndarray(float)`` or ``float``
            fitted y values
        """
        # ensure that self.opt_params has been set, indicating that fitting has already occurred
        if self.opt_params is None:
            msg = ('_CalibrationBase: transform: self.opt_params has not been set, self.fit must be successfully'
                   ' run prior to using this method')
            raise RuntimeError(msg)
        return self.fit_function(X, *self.opt_params)


class MassCalibration(_CalibrationBase):
    """
    mass calibration

    Attributes
    ----------
    mz_ref : ``numpy.ndarray(float)``
    mz_obs : ``numpy.ndarray(float)``
        arrays of reference and observed m/z values
    fit_func : ``str``
        function to use for fitting
    """
    # map valid fit functions to the actual functions and initial parameters
    _valid_fit_funcs = {
        'linear': (lambda X, a, b: a * X + b, (1., 0.)), 
    }

    def __init__(self, mz_ref, mz_obs, fit_func,
                 fit=True):
        """
        Initialize a new instance of MassCalibration

        Performs fitting at initialization.

        Parameters
        ----------
        mz_ref : ``numpy.ndarray(float)``
        mz_obs : ``numpy.ndarray(float)``
            arrays of reference and observed m/z values
        fit_func : ``str``
            function to use for fitting, valid options are: 'linear'
        fit : ``bool``, default=True
            perform fitting at initialization
        """
        # validate and store parameters
        self.mz_ref, self.mz_obs, self.fit_func  =  mz_ref, mz_obs, fit_func
        if self.fit_func not in self._valid_fit_funcs:
            msg = 'MassCalibration: __init__: fit_func "{}" invalid, must be one of: {}'
            valid_func_names = [_ for _ in self._valid_fit_funcs.keys()]
            raise ValueError(msg.format(self.fit_func, valid_func_names))
        # set fit_function and init_params based on fit_func
        self.fit_function, self.init_params = self._valid_fit_funcs[self.fit_func]
        if fit:
            # setup values for fitting, make corrections if applicable
            self._X = self.mz_obs
            self._y = self.mz_ref
            # fit calibration curve
            self._y_fit = self.fit(self._X, self._y)

    def calibrated_mass(self, mz):
        """
        returns calibrated masses for a set of uncalibrated masses

        Parameters
        ----------
        mz : ``numpy.ndarray(float)`` or ``float``
            uncalibrated m/z value(s)

        Returns
        -------
        mz_cal : ``numpy.ndarray(float)`` or ``float``
            calibrated m/z value(s)
        """
        return self.transform(mz)


class CCSCalibrationDTsf(_CalibrationBase):
    """
    single-field DT CCS calibration

    Attributes
    ----------
    mz : ``numpy.ndarray(float)``
        calibrant m/z values
    arrival_time : ``numpy.ndarray(float)``
        calibrant arrival times
    ref_ccs : ``numpy.ndarray(float)``
        calibrant reference CCS values
    z : ``float``
        charge state (converted to float)
    buffer_gas : ``str``
        buffer gas for IM separation
    """

    # map valid buffer gasses to their mass
    _valid_buffer_gasses = {
        'N2': monoiso_mass({'N': 2}), 
        'He': monoiso_mass({'He': 1}),
    }

    def __init__(self, mz, arrival_time, ref_ccs, z,
                 buffer_gas='N2', fit=True):
        """
        Initialize a new instance of CCSCalibrationDTsf 

        Performs fitting at initialization.

        Parameters
        ----------
        mz : ``numpy.ndarray(float)``
            calibrant m/z values
        arrival_time : ``numpy.ndarray(float)``
            calibrant arrival times
        ref_ccs : ``numpy.ndarray(float)``
            calibrant reference CCS values
        z : ``int``
            charge state
        buffer_gas : ``str``, default='N2'
            specify buffer gas for IM separation
        fit : ``bool``, default=True
            perform fitting at initialization
        """
        # validate and store parameters
        self.mz, self.arrival_time, self.ref_ccs, self.z =  mz, arrival_time, ref_ccs, float(z)
        # set init_params
        self.init_params = (0., 1.)  # t_fix = 0 to start, beta = 1 to start
        self.buffer_gas = buffer_gas
        if self.buffer_gas not in self._valid_buffer_gasses:
            msg = 'CCSCalibrationDTsf: __init__: buffer_gas "{}" invalid, must be one of: {}'
            raise ValueError(msg.format(self.buffer_gas, self._valid_buffer_gasses))
        self.buffer_gas_mass = self._valid_buffer_gasses[self.buffer_gas]
        if fit:
            # setup values for fitting
            self._X = np.array([self.arrival_time, self.mz])
            self._y = self.ref_ccs
            # fit calibration curve
            self._y_fit = self.fit(self._X, self._y)

    def fit_function(self, X, t_fix, beta):
        """
        X should have shape (2, ?) and contain vectors for arrival time and m/z
        function: CCS = z (arrival_time + t_fix) / (beta * m(mz)) 
        where t_fix and beta are the fit parameters
        and m(mz) is the mass term [sqrt(m_i / (m_i + m_b))] that accounts for ion and buffer gas masses
        """
        return self.z * (X[0] + t_fix) / (beta * self._mass_term(X[1]))

    def _mass_term(self, mz):
        """
        compute sqrt(m_i / (m_i + m_b)) term for m/z of ion and buffer gas
        """
        m = mz * self.z  # mass not m/z
        return np.sqrt(m / (m + self.buffer_gas_mass))

    def calibrated_ccs(self, mz, arrival_time):
        """
        returns calibrated CCS values for a set of m/z values and arrival times (also works for single values)

        Parameters
        ----------
        mz : ``numpy.ndarray(float)`` or ``float``
            m/z value(s)
        arrival_time : ``numpy.ndarray(float)`` or ``float``
            arrival time(s)

        Returns
        -------
        calibrated_ccs : ``numpy.ndarray(float)`` or ``float``
            calibrated CCS value(s)
        """
        if type(mz) == float:
            X = np.array([[arrival_time], [mz]])
        else:
            X = np.array([arrival_time, mz])
        return self.transform(X)


class CCSCalibrationTW(_CalibrationBase):
    """
    TWIMS CCS calibration

    Attributes
    ----------
    mz : ``numpy.ndarray(float)``
        calibrant m/z values
    arrival_time : ``numpy.ndarray(float)``
        calibrant arrival times
    ref_ccs : ``numpy.ndarray(float)``
        calibrant reference CCS values
    z : ``float``
        charge state (converted to float)
    fit_func : ``str``
        specify the type of function to use for fitting the calibration curve. Valid options are: ''
    correct_ccs : ``bool``
        perform CCS correction for charge state and reduced mass with buffer gas
    correct_dt : ``bool``
        perform arrival time correction for mass-dependent flight time outside of mobility region, this only
        really applies for cases where the mobility separation region is sufficiently small so that the time
        outside this region must be accounted for and is instrument-specific.
    buffer_gas : ``str``
        buffer gas for IM separation
    edc : ``float``
        EDC delay coefficient for arrival time correction (if used)
    """

    # map valid fit functions to the actual functions and initial parameters
    _valid_fit_funcs = {
        'linear': (lambda X, a, b: a * X + b, (1., 0.)), 
        'quadratic': (lambda X, a, b, c: a * X**2 + b * X + c, (0., 1., 0.)),
        'power1': (lambda X, a, b, c: a + b * np.power(X, c), (500., 1., 0.5)), 
        'power2': (lambda X, a, b, c: a * np.power((X + b), c), (1., 1e-4, 0.5)),
    }
    # map valid buffer gasses to their mass
    _valid_buffer_gasses = {
        'N2': monoiso_mass({'N': 2}), 
        'He': monoiso_mass({'He': 1}),
    }

    def __init__(self, mz, arrival_time, ref_ccs, z, fit_func,
                 correct_ccs=True, correct_dt=False, edc=None, buffer_gas='N2', fit=True):
        """
        Initialize a new instance of CCSCalibrationTW 

        Performs fitting at initialization.

        Parameters
        ----------
        mz : ``numpy.ndarray(float)``
            calibrant m/z values
        arrival_time : ``numpy.ndarray(float)``
            calibrant arrival times
        ref_ccs : ``numpy.ndarray(float)``
            calibrant reference CCS values
        z : ``int``
            charge state
        fit_func : ``str``
            specify the type of function to use for fitting the calibration curve. Valid options are: 'linear',
            'quadratic', 'power1', 'power2'
        correct_ccs : ``bool``, default=True
            perform CCS correction for charge state and reduced mass with buffer gas
        correct_dt : ``bool``, default=False
            perform arrival time correction for mass-dependent flight time outside of mobility region, this only
            really applies for cases where the mobility separation region is sufficiently small so that the time
            outside this region must be accounted for and is instrument-specific. If this is set to True, then the
            edc kwarg is also expected to be set.
        edc : ``float``, optional
            if arrival times should be corrected, set this as the EDC delay coefficient
        buffer_gas : ``str``, default='N2'
            specify buffer gas for IM separation
        fit : ``bool``, default=True
            perform fitting at initialization
        """
        # validate and store parameters
        self.mz, self.arrival_time, self.ref_ccs, self.z, self.fit_func =  mz, arrival_time, ref_ccs, float(z), fit_func
        if self.fit_func not in self._valid_fit_funcs:
            msg = 'CCSCalibrationTW: __init__: fit_func "{}" invalid, must be one of: {}'
            valid_func_names = [_ for _ in self._valid_fit_funcs.keys()]
            raise ValueError(msg.format(self.fit_func, valid_func_names))
        # set fit_function and init_params based on fit_func
        self.fit_function, self.init_params = self._valid_fit_funcs[self.fit_func]
        self.correct_ccs, self.correct_dt, self.edc = correct_ccs, correct_dt, edc
        if self.correct_dt and self.edc is None:
            msg = 'CCSCalibrationTW: __init__: correct_dt was set but no edc was provided'
            raise ValueError(msg)
        self.buffer_gas = buffer_gas
        if self.buffer_gas not in self._valid_buffer_gasses:
            msg = 'CCSCalibrationTW: __init__: buffer_gas "{}" invalid, must be one of: {}'
            raise ValueError(msg.format(self.buffer_gas, self._valid_buffer_gasses))
        self.buffer_gas_mass = self._valid_buffer_gasses[self.buffer_gas]
        if fit:
            # setup values for fitting, make corrections if applicable
            self._X = self._correct_dt(self.arrival_time, self.mz) if self.correct_dt else self.arrival_time
            self._y = self._correct_ccs(self.ref_ccs, self.mz) if self.correct_ccs else self.ref_ccs
            # fit calibration curve
            self._y_fit = self.fit(self._X, self._y)

    def _reduced_mass(self, mz):
        """
        compute reduced mass for an m/z with the buffer gas
        """
        m = mz * self.z  # mass not m/z
        return m * self.buffer_gas_mass / (m + self.buffer_gas_mass)
    
    def _correct_ccs(self, ccs, mz):
        """
        CCS' = CCS / (z * sqrt(1 / reduced_mass))
        """
        return ccs / (self.z * np.sqrt(1. / self._reduced_mass(mz)))

    def _inverse_correct_ccs(self, ccs_corr, mz):
        """
        CCS = CCS' * z * sqrt(1 / reduced_mass)
        """
        return ccs_corr * self.z * np.sqrt(1. / self._reduced_mass(mz))

    def _correct_dt(self, arrival_time, mz):
        """
        dt' = dt - (edc * sqrt(mz) / 1000)
        """
        return arrival_time - (self.edc * np.sqrt(mz) / 1e3)

    def calibrated_ccs(self, mz, arrival_time):
        """
        returns calibrated CCS values for a set of m/z values and arrival times (also works for single values)

        Parameters
        ----------
        mz : ``numpy.ndarray(float)`` or ``float``
            m/z value(s)
        arrival_time : ``numpy.ndarray(float)`` or ``float``
            arrival time(s)

        Returns
        -------
        calibrated_ccs : ``numpy.ndarray(float)`` or ``float``
            calibrated CCS value(s)
        """
        ccs = self.transform(self._correct_dt(arrival_time, mz) if self.correct_dt else arrival_time)
        return self._inverse_correct_ccs(ccs, mz) if self.correct_ccs else ccs
