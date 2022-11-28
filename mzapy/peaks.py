"""
mzapy/peaks.py

Dylan Ross (dylan.ross@pnnl.gov)

    module for peak fitting / feature finding / signal processing utilities
"""


import numpy as np
from scipy import signal, optimize, interpolate


def lerp_1d(x, y, x_min, x_max, density, threshold_y=True):
    """
    performs linear interpolation of arbitrary x, y data returning new x and y values
    between x_min and x_max at a specified density (# points / units of x)
    any interpolated values < 0 (e.g., when extrapolating beyond range of input values) get set to 0

    Parameters
    ----------
    x : ``list(float)``
        input x values
    y : ``list(float)``
        input y values
    x_min : ``float``
        minimum x value for interpolation
    x_max : ``float``
        maximum x value for interpolation
    density : ``float``
        interpolation density (# points / units of x)
    threshold_y : ``bool``, default=True
        ensure that no y values go below 0
    
    Returns
    -------
    x_interp : ``numpy.ndarray(float)``
        interpolated x values
    y_interp : ``numpy.ndarray(float)``
        interpolated y values
    """
    # x and y must be same length, and greater than 2
    if len(x) != len(y):
        msg = 'lerp_1d: x and y arrays must be same length'
        raise ValueError(msg)
    if len(x) < 2 or len(y) < 2:
        msg = 'lerp_1d: x and y arrays must have >= 2 points for interpolation'
        raise ValueError(msg)
    n_points = int(density * (x_max - x_min))
    if n_points < 1:
        msg = 'lerp_1d: x axis bounds ({}, {}) too small for specified density ({})'
        raise ValueError(msg.format(x_min, x_max, density))
    x_interp = np.linspace(x_min, x_max, n_points)
    y_interp = interpolate.interp1d(x, y, kind='linear', fill_value='extrapolate')(x_interp)
    # make sure no values below 0
    if threshold_y:
        y_interp[y_interp < 0] = 0
    return x_interp, y_interp


def lerp_2d(x, y, z, x_min, x_max, x_density, y_min, y_max, y_density, threshold_z=True):
    """
    performs linear interpolation of arbitrary x, y, z data returning new gridded x, y and z values
    between x_min, x_max, y_min, y_max at a specified density (# points / units of x, y)
    any interpolated values < 0 (e.g., when extrapolating beyond range of input values) get set to 0

    Parameters
    ----------
    x : ``list(float)``
        input x values
    y : ``list(float)``
        input y values
    z : ``list(float)``
        input z values
    x_min : ``float``
        minimum x value for interpolation
    x_max : ``float``
        maximum x value for interpolation
    x_density : ``float``
        x interpolation density (# points / units of x)
    y_min : ``float``
        minimum y value for interpolation
    y_max : ``float``
        maximum y value for interpolation
    y_density : ``float``
        y interpolation density (# points / units of y)
    threshold_z : ``bool``, default=True
        ensure that no z values go below 0
    
    Returns
    -------
    x_interp : ``numpy.ndarray(float)``
        interpolated x values
    y_interp : ``numpy.ndarray(float)``
        interpolated y values
    zz_interp : ``numpy.ndarray(float)``
        interpolated z values, 2D grid data 
    """
    n_points_x = int(x_density * (x_max - x_min))
    if n_points_x < 1:
        msg = '_lerp_2d: x axis bounds ({}, {}) too small for specified density ({})'
        raise ValueError(msg.format(x_min, x_max, x_density))
    n_points_y = int(y_density * (y_max - y_min))
    if n_points_y < 1:
        msg = '_lerp_2d: y axis bounds ({}, {}) too small for specified density ({})'
        raise ValueError(msg.format(y_min, y_max, y_density))
    x_interp = np.linspace(x_min, x_max, n_points_x)
    y_interp = np.linspace(y_min, y_max, n_points_y)
    xx_interp, yy_interp = np.meshgrid(x_interp, y_interp)
    z_interp = interpolate.griddata((x, y), z, (xx_interp, yy_interp), fill_value=0, method='linear')
    # make sure no z values below 0
    if threshold_z:
        z_interp[z_interp < 0] = 0
    return xx_interp, yy_interp, z_interp


def find_peaks_1d_localmax(x, y, min_rel_height, min_abs_height, fwhm_min, fwhm_max, min_dist):
    """
    find peaks in x, y data using local maximum method (scipy.signal.find_peaks)
    requires dense x data, uniformly spaced, monotonically increasing

    Parameters
    ----------
    x : ``np.array(float)``
        x data
    y : ``np.array(float)``
        y data
    min_rel_height : ``float``
        minimum height of peaks (relative to max y)
    min_abs_height : ``float``
        minimum absolute height
    fwhm_min : ``float``
        minimum peak width (FWHM)
    fwhm_max : ``float``
        maximum peak width (FWHM)
    min_dist : ``float``
        minimum distance (in units of x) between consecutive peaks

    Returns
    -------
    peak_means : ``numpy.ndarray(float)``
        mean x values for peaks 
    peak_heights : ``numpy.ndarray(float)``
        heights for peaks
    peak_fwhms : ``numpy.ndarray(float)``
        widths (FWHM) for peaks
    """
    # convert minimum distance (in units of x) to bins
    # compute average bin width
    bin_width = (x[-1] - x[0]) / len(x)
    min_bins = int(min_dist / bin_width)
    min_bins = 1 if min_bins < 1 else min_bins
    fwhm_min_bins, fwhm_max_bins = int(fwhm_min / bin_width), int(fwhm_max / bin_width)
    x, y = np.array(x), np.array(y)
    min_height = max(min_rel_height * np.max(y), min_abs_height)
    peak_idxs, peak_props = signal.find_peaks(y, height=min_height, distance=min_bins, 
                                              width=(fwhm_min_bins, fwhm_max_bins))
    peak_means = x[peak_idxs]
    peak_heights = peak_props['peak_heights']
    peak_fwhms = bin_width * signal.peak_widths(y, peak_idxs, rel_height=0.5)[0]
    return peak_means, peak_heights, peak_fwhms


def _gauss(x, mean, height, fwhm):
    """
    gaussian function for peak fitting
    FWHM converted to standard deviation
    y = height * exp(-(x - mean)^2 / (0.3606 * FWHM^2))

    Parameters
    ----------
    x : ``float`` or ``numpy.ndarray(float)``
        input data 
    mean : ``float``
        mean of gaussian distribution
    height : ``float``
        amplitude of gaussian distribution
    fwhm : ``float``
        FWHM of gaussian distribution (derived from standard deviation)

    Returns
    -------
    y : ``float`` or ``numpy.ndarray(float)``
        output data
    """
    return height * np.exp(-(x - mean)**2. / (0.3606 * fwhm**2.))


# returns peak params (mean, height, fwhm) and residual y (values < 0 set to 0)
def _gauss_fit_one_peak(x, y, fwhm_min, fwhm_max, truncate_y_resids):
    """
    performs gaussian fitting for a single peak from x, y data 

    Parameters
    ----------
    x : ``numpy.ndarray(float)``
        x data
    y : ``numpy.ndarray(float)``
        y data
    fwhm_min : ``float``
        minimum peak width (FWHM)
    fwhm_max : ``float``
        maximum peak width (FWHM)
    truncate_y_resids : ``bool``
        ensure that residuals do not have any intensity values below 0

    Returns
    -------
    params : ``tuple(float)``
        fitted gaussian function parameters (mean, height, fwhm) defining the peak
    y_resids : ``numpy.ndarray(float)``
        fit residuals for y data
    """
    # initial guess of parameters
    mean_0 = x[np.argmax(y)]
    height_0 = np.max(y)
    fwhm_0 = fwhm_min * 1.05 # initial FWHM is min FWHM * 1.05
    bounds = ([np.min(x), 0, fwhm_min], [np.max(x), np.inf, fwhm_max])
    try:
        params, _ = optimize.curve_fit(_gauss, x, y, p0=(mean_0, height_0, fwhm_0), bounds=bounds, maxfev=5000)
    except RuntimeError:
        # least-squares fit failed
        return None, None
    y_resid = y - _gauss(x, *params)
    if truncate_y_resids:
        # make sure no values below 0
        y_resid[y_resid < 0] = 0
    return params, y_resid


def find_peaks_1d_gauss(x, y, min_rel_height, min_abs_height, fwhm_min, fwhm_max, max_peaks, truncate_y_resids):
    """
    find peaks in x, y data using a successive gaussian function fit method
    requires dense x data, uniformly spaced, monotonically increasing. the most intense 
    peak in the data is fitted using a gaussian function, then the residuals are fitted for
    the next peak. This process is repeated until a stopping criterion is reached:

    * peak height is lower than a specified absolute or relative threshold
    * a maximum number of peaks have been found

    Parameters
    ----------
    x : ``numpy.ndarray(float)``
        x data
    y : ``numpy.ndarray(float)``
        y data
    min_rel_height : ``float``
        minimum height of peaks (relative to max y)
    min_abs_height : ``float``
        minimum absolute height
    fwhm_min : ``float``
        minimum peak width (FWHM)
    fwhm_max : ``float``
        maximum peak width (FWHM)
    max_peaks : ``int``
        maximum number of peaks to find
    truncate_y_resids : ``bool``
        ensure that residuals do not have any intensity values below 0

    Returns
    -------
    peak_means : ``numpy.ndarray(float)``
        mean x values for peaks 
    peak_heights : ``numpy.ndarray(float)``
        heights for peaks
    peak_fwhms : ``numpy.ndarray(float)``
        widths (FWHM) for peaks
    """
    peak_means, peak_heights, peak_fwhms = [], [], []
    min_height = max(min_abs_height, min_rel_height * max(y))
    n = 0
    while True:
        # try fitting a peak
        peak, peak_resids = _gauss_fit_one_peak(x, y, fwhm_min, fwhm_max, truncate_y_resids)
        if peak is not None:
            if peak[1] >= min_height and n < max_peaks:
                # peak fit is acceptable
                peak_means.append(peak[0])
                peak_heights.append(peak[1])
                peak_fwhms.append(peak[2])
                n += 1
                y = peak_resids
            else:
                # peak is not intense enough or n >= max_peaks
                break
        else:
            # peak fitting failed
            break
    return np.array(peak_means), np.array(peak_heights), np.array(peak_fwhms)


def calc_gauss_psnr(x, y, peak_params):
    """
    Computes an estimate of the peak signal to noise ratio (pSNR) by comparing the peak height of a fitted peak
    (valid only for gaussian peak fitting method) to the RMS of the residuals from the peak fitting. *This method of 
    estimating the pSNR is sensitive to the intensity of the peak relative to background noise AND the goodness of fit 
    of the gaussian function. Any deviation of the signal from the ideal gaussian shape and/or presence of multiple 
    peaks will decrease the pSNR, even if the peak height is much greater than the true magnitude of noise.*

    Parameters
    ----------
    x : ``numpy.ndarray(float)``
        x data
    y : ``numpy.ndarray(float)``
        y data
    peak_params : ``tuple(float)``
        gaussian function fit parameters (mean, height, fwhm) defining the peak

    Returns
    -------
    psnr : ``float``
        estimate of the peak signal to noise ratio
    """
    pdt, pht, pwt = peak_params
    resids = _gauss(x, pdt, pht, pwt) - y
    psnr = pht / np.sqrt(np.mean(resids * resids))
    return psnr


def calc_peak_area(peak_height, peak_fwhm):
    """
    compute the peak area from peak height and FWHM, assumes gaussian peak shape (also works with values from localmax 
    method, but may be less accurate)

    Parameters
    ----------
    peak_height : ``float``
        peak height
    peak_fwhm : ``float``
        peak FWHM

    Returns
    -------
    area : ``float``
        peak area
    """
    # area = sqrt(pi / ln(2)) * h * w / 2  ~= 1.064467 * h * w
    return 1.064467 * peak_height * peak_fwhm

