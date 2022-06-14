"""
mzapy/calibration.py
    
Dylan Ross (dylan.ross@pnnl.gov)

    Module for performing calibration operations.

    CCS for TWIM measurements, method based on:
        - Anal. Chem. 2016, 88, 7329−7336 (calibration power function)
        - Nat. Protoc. 2008, 3, 1139–1152 (dt and CCS correction)
"""


import math

from scipy import optimize
from matplotlib import pyplot as plt, gridspec as gs
import numpy as np

from mzapy.isotopes import monoiso_mass


class _CCSCalibrationBase:
    """
    Modular base object for making and applying CCS calibrations, not meant to be used directly

    Attributes
    ----------
    ref_mass : ``float``
        reference mass (for reduced mass calculation)
    edc : ``float``
        EDC delay coefficient
    charge : ``float``
        charge state
    calibrated : ``bool``
        indicates whether a calibration has set up
    cal_params : ``tuple(float or None)``
        calibration parameters (A, t0, B) or (None, None, None) if calibration has not been set up
    cal_data : ``tuple(list(float) or None)``
        calibration data (mz, dt, ref_ccs) or (None, None, None) if calibration has not been set up
    cal_corr_data : ``tuple(list(float) or None)``
        corrected calibration data (corr_dt, corr_ref_ccs) or (None, None) if calibration has not been set up
    """

    def __init__(self, ref_mass, edc, charge):
        """
        inits a new instance of CCSCalibrationBase

        Parameters
        ----------
        ref_mass : ``float``
            reference mass (drift gas)
        edc : ``float``
            EDC delay coefficient
        charge : ``int``
            charge state
        """
        self.ref_mass = ref_mass
        self.edc = edc
        self.charge = float(charge)
        # flag whether a calibration curve has been fit
        self.calibrated = False
        # calibration curve parameters
        self.cal_params = (None, None, None)
        # masses, dts, and reference ccss used to fit the calibration curve
        self.cal_data = (None, None, None)
        # corrected dts and reference ccss
        self.cal_corr_data = (None, None)

    def _reduced_mass(self, mz):
        """
        Calculates reduced mass of an ion using a reference mass (self.ref_mass)
        
        Paramters
        ---------
        mz : ``float``
            input m/z
        
        Returns
        -------
        reduced_mass : ``float``
            reduced mass
        """
        return (mz * self.ref_mass) / (mz + self.ref_mass)

    def _corrected_dt(self, dt, mz):
        """
        Calculates a drift time corrected for mass-dependent flight time outside of mobility region
        uses EDC delay coefficient in self.edc, 
        alternative strategy is to calculate this flight time using transfer region length and wave velocity
        
        Paramters
        ---------
        dt : ``float``
            original drift time
        mz : ``float``
            m/z to use for the correction

        Returns
        -------
        corr_dt : ``float``
            corrected drift time
        """
        return dt - (math.sqrt(mz) * self.edc) / 1000.

    def _corrected_ccs(self, ccs, mz):
        """
        Calculates CCS corrected for mass-dependent flight time

        Parameters
        ----------
        ccs : ``float``
            original CCS
        mz : ``float``
            m/z to use for the correction
        
        Returns
        -------
        corr_ccs : ``float``
            corrected CCS
        """
        return ccs * (self.charge * math.sqrt(self._reduced_mass(mz)))

    def _cal_curve(self, corr_dt, A, t0, B):
        """
        Basic power function for calibration curve of the form:
            CCS' = z * A * (dt' + t0)^B
        
        Parameters
        ----------
        corr_dt : ``float``
            corrected drift time
        A : ``float``
            A parameter
        t0 : ``float``
            t0 parameter
        B : ``float``
            B parameter

        Returns
        -------
        corr_ccs : ``float``
            corrected CCS
    """
        return self.charge * A * (corr_dt + t0)**B

    def _fit_cal_curve(self, mzs, dts, ccss):
        """
        Fits a calibration curve to dt and CCS values and stores the optimized
        parameters in ``self.cal_params``
        *Expects UNCORRECTED dts and ccss*

        Parameters
        ----------
        mzs : ``list(float)``
            calibrant m/zs
        dts : ``list(float)``
            uncorrected calibrant drift times
        ccss : ``list(float)``
            uncorrected calibrant ccss
        """
        # need to create a temporary instance of this class so self.cal_curve can be used as a bound method in curve_fit
        temp_cal = _CCSCalibrationBase(self.ref_mass, self.edc, self.charge)
        corr_dt = [self._corrected_dt(dt, mass) for dt, mass in zip(dts, mzs)]
        corr_ccs = [self._corrected_ccs(ccs, mass) for ccs, mass in zip(ccss, mzs)]
        self.cal_params, cov = optimize.curve_fit(temp_cal.cal_curve,
                                                  corr_dt, corr_ccs,
                                                  maxfev=1000000, p0=(500., 0.001, 0.5))
        # after a successful fit store the masses, drift times, ccss used to
        # create the calibration curve
        self.calibrated = True
        self.cal_data = (mzs, dts, ccss)
        self.cal_corr_data = (corr_dt, corr_ccs)

    def calibrated_ccs(self, mz, dt):
        """
        Use the calibration curve parameters in ``self.cal_params`` to
        compute calibrated CCS for an m/z dt pair
        
        Parameters
        ----------
        mz : ``float``
            m/z
        dt : ``float``
            drift time
        
        Returns
        -------
        ccs : ``float``
            calibrated ccs
        """
        if not self.calibrated:
            raise RuntimeError("CCSCalibrationBase: calibrated_ccs: calibration curve has not been fit yet")
        return self._cal_curve(self._corrected_dt(dt, mz), *self.cal_params) / (self.charge * math.sqrt(self._reduced_mass(mz)))

    def cal_curve_figure(self, fig_name):
        """
        Produces a figure from the CCS calibration curve with residuals and saves it as a .png
        
        Parameters
        ----------
        fig_name : ``str``
            filename to save the calibration curve under
        """
        if not self.calibrated:
            raise RuntimeError("CCSCalibrationBase: cal_curve_figure: calibration curve has not been fit yet")
        # gather necessary data
        corr_dt, corr_ccs = self.cal_corr_data
        masses, dts, ccss = self.cal_data
        corr_ccs_calc = [self._corrected_ccs(self.calibrated_ccs(mass, dt), mass) for mass, dt in zip(masses, dts)]
        ccs_resid = [100. * (ccs - self.calibrated_ccs(mass, dt)) / ccs for mass, dt, ccs in zip(masses, dts, ccss)]
        # set up the figure
        fig = plt.figure(figsize=(5, 4.5))
        grid = gs.GridSpec(2, 1, height_ratios=(5, 2))
        ax1, ax2 = fig.add_subplot(grid[0]), fig.add_subplot(grid[1])
        ax2.axhline(lw=0.75, ls="--", c="k")
        # plot the data
        ax1.plot(corr_dt, corr_ccs, "bo", ms=4, mew=1, fillstyle='none', label='calibrants')
        x = np.linspace(min(corr_dt), max(corr_dt), 100)
        y = self._cal_curve(x, *self.cal_params)
        ax1.plot(x, y, "b-", lw=1, alpha=0.6, label='fitted')
        ax1.legend(frameon=False)
        ax2.bar(corr_dt, ccs_resid, 0.25, color=(0., 0., 1., 0.5), align='center')
        # axis adjustments and labels
        ax1.set_title("CCS calibration")
        ax1.set_ylabel("corrected CCS")
        ax2.set_ylabel("residual CCS (%)")
        ax2.set_xlabel("corrected drift time (ms)")
        for d in ['top', 'right']:
            ax1.spines[d].set_visible(False)
            ax2.spines[d].set_visible(False)
        ax2.spines['bottom'].set_visible(False)
        plt.savefig(fig_name, dpi=400, bbox_inches='tight')
        plt.close()

    def __str__(self):
        """
        produces a string representation of this instance, contains information
        about the calibration curve (if fitted)
        
        Returns
        -------
        s : ``str``
            string representation of this CCSCalibrationBase object
        """
        if self.cal_params[0] is None:
            return "CCS calibration curve not fitted"
        # calibration has been completed, include the fitted calibration curve parameters
        else:
            return "CCS calibration fitted parameters: \n\tA = {:.3f}\n\tt0 = {:.3f}\n\tB = {:.3f}".format(*self.cal_params)


class CCSCalibrationList(_CCSCalibrationBase):
    """
    object for making and applying CCS calibrations using lists of m/z, dt, and reference CCS

    Attributes
    ----------
    ref_mass : ``float``
        reference mass (for reduced mass calculation)
    edc : ``float``
        EDC delay coefficient
    charge : ``float``
        charge state
    calibrated : ``bool``
        indicates whether a calibration has set up
    cal_params : ``tuple(float or None)``
        calibration parameters (A, t0, B) or (None, None, None) if calibration has not been set up
    cal_data : ``tuple(list(float) or None)``
        calibration data (mz, dt, ref_ccs) or (None, None, None) if calibration has not been set up
    cal_corr_data : ``tuple(list(float) or None)``
        corrected calibration data (corr_dt, corr_ref_ccs) or (None, None) if calibration has not been set up
    """

    def __init__(self, edc, charge, mz, dt, ref_ccs):
        """
        inits a CCS calibration object from lists of m/z, drift times, and reference CCS values for calibrants
        Automatically tries to fit calibration curve 
        Assumes N2 as drift gas for reduced mass calculation
        
        Parameters
        ----------
        edc : ``float``
            EDC delay coefficient
        charge : ``int``
            charge state
        mz : list(float)
            calibrant m/z values
        dt : ``list(float)``
            calibrant drift times
        ref_ccs : ``list(float)``
            calibrant reference CCS values
        """
        # run superclass __init__ with ref mass, EDC, and charge state
        # default to nitrogen for reduced mass calculation
        super().__init__(monoiso_mass({'N': 2}), edc, charge)

        # fit calibration curve
        self._fit_cal_curve(mz, dt, ref_ccs)

