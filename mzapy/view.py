"""
mzapy/view.py

Dylan Ross (dylan.ross@pnnl.gov)

    module for interactive viewing and/or generating plots of extracted data

    TODO (Dylan Ross): add some functions for plotting 2D data like RTMZ, DTMZ, RTDT profiles
"""


from functools import wraps

from matplotlib import pyplot as plt, rcParams

from mzapy._util import _ppm_error


# always set font size to 6 on import
rcParams['font.size'] = 6


def _setup_and_save_or_show_plot(label, ax, figname, figsize, ign_newax):
    """
    Decorator with boilerplate for setting up then saving or showing the resulting plot. Wraps a function that plots
    data and applies formatting. Wrapped functions should always return the ``Axes`` instance. The parameters of 
    the decorator set the default values for these kwargs.

    .. code-block:: python3
        :caption: template for ``plot_and_format`` function

        @_setup_and_save_or_show_plot(label=None, ax=None, figname=None, figsize=(3.33, 1.5), ign_newax=False)
        def plot_spectrum(mz, i, 
                          mz_range=None, c='b', 
                          **kwargs):
            # ----------
            # Parameters
            # ----------
            # xdata : ``numpy.ndarray(float)``
            # ydata : ``numpy.ndarray(float)``
            #     x and y data    
            #  ... other params ...
            # label : ``str``, optional
            #     if provided, label for the data being plotted
            # ax : ``matplotlib.axes.Axes``, optional
            #     if provided, do not make new ``Figure`` and ``Axes`` instances, just plot the spectrum and format it 
            #     within the provided ``Axes`` instance
            # figname : ``str``, optional
            #     if not provided, do not save the figure or display it for interactive viewing. if figname is "show" 
            #     then display the figure for interactive viewing, and if figname is a path to an image file then save 
            #     the image to that path
            # figsize : ``tuple(float)``, default=(3.33, 1.5)
            #     figure width, height in inches. Only used when ``ax`` is not set
            # ign_newax : ``bool``, default=False
            #     when True, ignores the newax flag and performs all plot customizations regardless of whether the existing
            #     ``matplotlib.axes.Axes`` instance has already been customized
            # -------
            # Returns
            # -------
            # ax : ``matplotlib.axes.Axes`` or ``None``
            #     ``Axes`` instance for adding further customizations to the plot. Returns None if figname was 
            #     provided since ``plt.close()`` gets called after interactive viewing or saving the plot so the 
            #     ``Axes`` instance cannot be used after that
            # -------

            # flags set by decorator
            newax = kwargs['_newax']
            ign_newax = kwargs['_ign_newax']
            ax = kwargs['ax']
            label = kwargs['label']
            ax.plot(xdata, ydata, ...)
            # apply formatting to the plot if the Axes instance is new, otherwise assume it has already been done
            # unless ign_newax is True, then ignore the newax flag
            if newax or ign_newax:
                ax.set_xlabel('X data')
                ax.set_ylabel('Y data')
                # take off the top and right borders
                for d in ['top', 'right']:
                    ax.spines[d].set_visible(False)
            # always return the Axes instance, decorator will decide whether to return Axes instance or None
            return ax


    Parameters
    ----------
    label : ``str``, optional
        if provided, label for the data being plotted
    ax : ``matplotlib.axes.Axes``, optional
        if provided, do not make new ``Figure`` and ``Axes`` instances, just plot the data and format it within the
        provided ``Axes`` instance
    figname : ``str``, optional
        if not provided, do not save the figure or display it for interactive viewing. if figname is "show" then
        display the figure for interactive viewing, and if figname is a path to an image file then save the image to
        that path
    figsize : ``tuple(float)``, default=(3.33, 2.22)
        figure width, height in inches. Only used when ``ax`` is not set
    ign_newax : ``bool``, default=False
        when True, ignores the newax flag and performs all plot customizations regardless of whether the existing
        ``matplotlib.axes.Axes`` instance has already been customized
    """
    def wrapper_outer(plot_and_format):
        # outer wrapper to handle kwargs with default values
        @wraps(plot_and_format)
        def wrapper(*args, label=label, ax=ax, figname=figname, figsize=figsize, ign_newax=ign_newax, **kwargs):
            # pre-plotting setup
            # flag indicating that a new Axes instance was made
            newax = False
            if ax is None:
                fig, ax = plt.subplots(figsize=figsize)
                newax = True
            # pass along the new axes instance, newax flag, and label
            kwargs['ax'] = ax
            kwargs['_newax'] = newax
            kwargs['_ign_newax'] = ign_newax
            kwargs['label'] = label
            # ------------------------------------------- 
            ax = plot_and_format(*args, **kwargs)
            # -------------------------------------------
            # save/show the plot if figname was provided, otherwise return Axes instance
            if figname is not None:
                # at this point all the data has been added to the plot that will be added, so now we can add the legend 
                # (if a label was provided)
                if label is not None:
                    ax.legend(frameon=False)
                plt.tight_layout()
                if figname == 'show':
                    plt.show()
                else:
                    plt.savefig(figname, dpi=350, bbox_inches='tight')
                plt.close()
                return None  # the plot is closed, no need to return an Axes instance
            else:
                # return the Axes instance so further customization can be done to it
                return ax
        return wrapper
    return wrapper_outer


@_setup_and_save_or_show_plot(label=None, ax=None, figname=None, figsize=(3.33, 1.5), ign_newax=False)
def plot_spectrum(mz, i, 
                  mz_range=None, c='b', ls='-',
                  **kwargs):
    """
    generate a plot of a mass spectrum (mz vs. i)

    Parameters
    ----------
    mz : ``numpy.ndarray(float)``
    i : ``numpy.ndarray(float)``
        m/z and intensity components of the mass spectrum
    mz_range : tuple(float), optional
        if provided, sets m/z bounds to display from spectrum
    c : ``str``, default='b'
        color to use for plotting spectrum, default is blue
    ls : ``str``, default='-'
        linestyle for plot, default is solid line
    label : ``str``, optional
        if provided, label for the spectrum being plotted
    ax : ``matplotlib.axes.Axes``, optional
        if provided, do not make new ``Figure`` and ``Axes`` instances, just plot the spectrum and format it within the
        provided ``Axes`` instance
    figname : ``str``, optional
        if not provided, do not save the figure or display it for interactive viewing. if figname is "show" then
        display the figure for interactive viewing, and if figname is a path to an image file then save the image to
        that path
    figsize : ``tuple(float)``, default=(3.33, 1.5)
        figure width, height in inches. Only used when ``ax`` is not set
    ign_newax : ``bool``, default=False
        when True, ignores the newax flag and performs all plot customizations regardless of whether the existing
        ``matplotlib.axes.Axes`` instance has already been customized

    Returns
    -------
    ax : ``matplotlib.axes.Axes`` or ``None``
        ``Axes`` instance for adding further customizations to the plot. Returns None if figname was provided since
        ``plt.close()`` gets called after interactive viewing or saving the plot so the ``Axes`` instance cannot be 
        used after that
    """
    # flags set by decorator
    newax = kwargs['_newax']
    ign_newax = kwargs['_ign_newax']
    ax = kwargs['ax']
    label = kwargs['label']
    # plot the spectrum
    ax.plot(mz, i, ls=ls, c=c, lw=0.75, label=label)
    # apply formatting to the plot if the Axes instance is new, otherwise assume it has already been done
    # unless ign_newax is True, then ignore the newax flag
    if newax or ign_newax:
        if mz_range is not None:
            ax.set_xlim(mz_range)
        ax.set_xlabel('m/z')
        ax.set_ylabel('intensity')
        for d in ['top', 'right']:
            ax.spines[d].set_visible(False)
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    # always return the Axes instance, decorator will decide whether to return Axes instance or None
    return ax


@_setup_and_save_or_show_plot(label=None, ax=None, figname=None, figsize=(3.33, 1.5), ign_newax=False)
def plot_chrom(rt, i, 
               rt_range=None, c='b', ls='-', rt_unit='min',
               **kwargs):
    """
    generate a plot of a chromatogram (rt vs. i)

    Parameters
    ----------
    rt : ``numpy.ndarray(float)``
    i : ``numpy.ndarray(float)``
        retention time and intensity components of the chromatogram
    rt_range : tuple(float), optional
        if provided, sets retention time bounds to display from chromatogram
    c : ``str``, default='b'
        color to use for plotting chromatogram, default is blue
    ls : ``str``, default='-'
        linestyle for plot, default is solid line
    rt_unit : ``str``, default='min'
        units of retention time, default is minutes (min)
    label : ``str``, optional
        if provided, label for the chromatogram being plotted
    ax : ``matplotlib.axes.Axes``, optional
        if provided, do not make new ``Figure`` and ``Axes`` instances, just plot the spectrum and format it within the
        provided ``Axes`` instance
    figname : ``str``, optional
        if not provided, do not save the figure or display it for interactive viewing. if figname is "show" then
        display the figure for interactive viewing, and if figname is a path to an image file then save the image to
        that path
    figsize : ``tuple(float)``, default=(3.33, 1.5)
        figure width, height in inches. Only used when ``ax`` is not set
    ign_newax : ``bool``, default=False
        when True, ignores the newax flag and performs all plot customizations regardless of whether the existing
        ``matplotlib.axes.Axes`` instance has already been customized

    Returns
    -------
    ax : ``matplotlib.axes.Axes`` or ``None``
        ``Axes`` instance for adding further customizations to the plot. Returns None if figname was provided since
        ``plt.close()`` gets called after interactive viewing or saving the plot so the ``Axes`` instance cannot be 
        used after that
    """
    # flags set by decorator
    newax = kwargs['_newax']
    ign_newax = kwargs['_ign_newax']
    ax = kwargs['ax']
    label = kwargs['label']
    # plot the spectrum
    ax.plot(rt, i, ls=ls, c=c, lw=0.75, label=label)
    # apply formatting to the plot if the Axes instance is new, otherwise assume it has already been done
    # unless ign_newax is True, then ignore the newax flag
    if newax or ign_newax:
        if rt_range is not None:
            ax.set_xlim(rt_range)
        ax.set_xlabel('retention time ({})'.format(rt_unit))
        ax.set_ylabel('intensity')
        for d in ['top', 'right']:
            ax.spines[d].set_visible(False)
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    # always return the Axes instance, decorator will decide whether to return Axes instance or None
    return ax


@_setup_and_save_or_show_plot(label=None, ax=None, figname=None, figsize=(3.33, 1.5), ign_newax=False)
def plot_atd(dt, i, 
             dt_range=None, c='b', ls='-', dt_unit='ms',
             **kwargs):
    """
    generate a plot of an arrival time distribution (dt vs. i)

    Parameters
    ----------
    dt : ``numpy.ndarray(float)``
    i : ``numpy.ndarray(float)``
        drift time and intensity components of the arrival time distribution
    dt_range : tuple(float), optional
        if provided, sets drift time bounds to display from arrival time distribution
    c : ``str``, default='b'
        color to use for plotting arrival time distribution, default is blue
    ls : ``str``, default='-'
        linestyle for plot, default is solid line
    dt_unit : ``str``, default='ms'
        units of drift time, default is milliseconds (ms)
    label : ``str``, optional
        if provided, label for the arrival time distribution being plotted
    ax : ``matplotlib.axes.Axes``, optional
        if provided, do not make new ``Figure`` and ``Axes`` instances, just plot the spectrum and format it within the
        provided ``Axes`` instance
    figname : ``str``, optional
        if not provided, do not save the figure or display it for interactive viewing. if figname is "show" then
        display the figure for interactive viewing, and if figname is a path to an image file then save the image to
        that path
    figsize : ``tuple(float)``, default=(3.33, 1.5)
        figure width, height in inches. Only used when ``ax`` is not set
    ign_newax : ``bool``, default=False
        when True, ignores the newax flag and performs all plot customizations regardless of whether the existing
        ``matplotlib.axes.Axes`` instance has already been customized

    Returns
    -------
    ax : ``matplotlib.axes.Axes`` or ``None``
        ``Axes`` instance for adding further customizations to the plot. Returns None if figname was provided since
        ``plt.close()`` gets called after interactive viewing or saving the plot so the ``Axes`` instance cannot be 
        used after that
    """
    # flags set by decorator
    newax = kwargs['_newax']
    ign_newax = kwargs['_ign_newax']
    ax = kwargs['ax']
    label = kwargs['label']
    # plot the spectrum
    ax.plot(dt, i, ls=ls, c=c, lw=0.75, label=label)
    # apply formatting to the plot if the Axes instance is new, otherwise assume it has already been done
    # unless ign_newax is True, then ignore the newax flag
    if newax or ign_newax:
        if dt_range is not None:
            ax.set_xlim(dt_range)
        ax.set_xlabel('arrival time ({})'.format(dt_unit))
        ax.set_ylabel('intensity')
        for d in ['top', 'right']:
            ax.spines[d].set_visible(False)
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    # always return the Axes instance, decorator will decide whether to return Axes instance or None
    return ax


def add_peaks_to_plot(ax, peak_means, peak_heights, peak_fwhms, 
                      c='r', add_text_lbl=False, x_units='', fontsize=6):
    """
    annotate an existing plot with peaks showing their means, heights, and FWHMs as crosses, can optionally add a 
    text label with the format {mean} +/- {FWHM / 2} {x_units} ({height})
    
    Parameters
    ----------
    ax : ``matplotlib.axes.Axes``
        ``Axes`` instance containing plot to annotate with peaks
    peak_means : ``list(float)``
    peak_heights : ``list(float)``
    peak_fwhms : ``list(float)``
        peak parameters (mean, height, FWHM) as lists
    c : ``str``, default='r'
        line color for annotated peaks
    add_text_lbl : ``bool``, default=False
        whether to annotate peaks with a text label, format: {mean} +/- {FWHM / 2} {x_units} ({height})
    x_units : ``str``, default=''
        when peaks are annotated with text labels, add specified units to the part of the label with the mean and FWHM
    fontsize : ``int``, default=6
        when peaks are annotated with text labels, specify font size for text labels
    """
    for pk_mean, pk_height, pk_width in zip(peak_means, peak_heights, peak_fwhms):
        # annotate with a cross showing the peak mean, height, and FWHM
        ax.plot([pk_mean - (pk_width / 2), pk_mean + (pk_width / 2)], [pk_height / 2, pk_height / 2], 
                ls='-', c=c, lw=1)
        ax.plot([pk_mean, pk_mean], [0, pk_height], 
                ls='-', c=c, lw=1)
        if add_text_lbl:
            s = '{:.1e} +/- {:.1e} {}\n{:.1e}'.format(pk_mean, pk_width, x_units, pk_height)
            # add some text with the peak info
            ax.text(pk_mean + pk_width / 2, 0.75 * pk_height, s, size=fontsize, c=c)


def plot_mass_calibration(calibration, figname=None):
    """
    Generates a plot of a mass calibration from a fitted instance of ``mzapy.calibration.MassCalibration``

    Parameters
    ----------
    calibration : ``mzapy.calibration.MassCalibration``
        fitted mass calibration instance
    figname : ``str``, optional
        if not provided, do not save the figure or display it for interactive viewing. if figname is "show" then
        display the figure for interactive viewing, and if figname is a path to an image file then save the image to
        that path
    """
    # ensure that the calibration has been fitted
    if calibration.opt_params is None:
        msg = 'plot_mass_calibration: calibration has not been fit yet (calibration.opt_params is None)'
        raise RuntimeError(msg)
    # get reference and fitted values
    X, y, y_fit = calibration._X, calibration._y, calibration._y_fit
    ppms = _ppm_error(y, y_fit)
    # get a text description of the fit parameters
    text_format = {
        'linear': "Y = {:.6f} * X {:+.6f}",
    }
    param_text = text_format[calibration.fit_func].format(*calibration.opt_params)
    fig, axs = plt.subplots(nrows=2, figsize=(3.33, 4), gridspec_kw={'height_ratios': [2.5, 1]})
    xlab = 'm/z (obs.)'
    axs[0].plot(X, y, 'b.', ms=5)
    axs[0].plot(X, y_fit, 'b--', lw=0.75)
    axs[0].set_ylabel('m/z (ref.)')
    axs[0].text(0.05, 1., param_text, color='b', transform=axs[0].transAxes)
    axs[1].plot(X, ppms, 'b.', ms=5)
    axs[1].axhline(0, ls='--', lw=0.5, color='k')
    axs[1].set_xlabel(xlab)
    axs[1].set_ylabel('ppm error')
    for ax in axs:
        for d in ['top', 'right']:
            ax.spines[d].set_visible(False)
    if figname is not None:
        plt.tight_layout()
        if figname == 'show':
            plt.show()
        else:
            plt.savefig(figname, dpi=350, bbox_inches='tight')


def plot_dtsf_ccs_calibration(calibration, figname=None):
    """
    Generates a plot of a single-filed DTIMS CCS calibration from a fitted 
    instance of ``mzapy.calibration.CCSCalibrationDTsf``

    Parameters
    ----------
    calibration : ``mzapy.calibration.CCSCalibrationDTsf``
        fitted single-field DTIMS CCS calibration instance
    figname : ``str``, optional
        if not provided, do not save the figure or display it for interactive viewing. if figname is "show" then
        display the figure for interactive viewing, and if figname is a path to an image file then save the image to
        that path
    """
    # ensure that the calibration has been fitted
    if calibration.opt_params is None:
        msg = 'plot_dtsf_ccs_calibration: calibration has not been fit yet (calibration.opt_params is None)'
        raise RuntimeError(msg)
    # get reference and fitted values
    y, y_fit = calibration._y, calibration._y_fit
    dt = calibration.arrival_time
    ccs_ref = calibration.ref_ccs
    ccs_cal = calibration.calibrated_ccs(calibration.mz, calibration.arrival_time)
    percent_error = 100. * (ccs_cal - ccs_ref) / ccs_ref
    # get a text description of the fit parameters
    text_format = 't_fix = {:+.4f}\nbeta = {:.4f}' 
    param_text = text_format.format(*calibration.opt_params)
    fig, axs = plt.subplots(nrows=2, figsize=(3.33, 4), gridspec_kw={'height_ratios': [2.5, 1]})
    xlab = 'arrival time (ms)'
    axs[0].plot(dt, y, 'b.', ms=5)
    axs[0].plot(dt, y_fit, 'b--', lw=0.75)
    axs[0].set_ylabel('CCS (Ang^2)')
    axs[0].text(0.05, 1., param_text, color='b', transform=axs[0].transAxes)
    axs[1].plot(dt, percent_error, 'b.', ms=5)
    axs[1].axhline(0, ls='--', lw=0.5, color='k')
    axs[1].set_xlabel(xlab)
    axs[1].set_ylabel('CCS % error')
    for ax in axs:
        for d in ['top', 'right']:
            ax.spines[d].set_visible(False)
    if figname is not None:
        plt.tight_layout()
        if figname == 'show':
            plt.show()
        else:
            plt.savefig(figname, dpi=350, bbox_inches='tight')


def plot_tw_ccs_calibration(calibration, figname=None):
    """
    Generates a plot of a TWIMS CCS calibration from a fitted 
    instance of ``mzapy.calibration.CCSCalibrationTW``

    Parameters
    ----------
    calibration : ``mzapy.calibration.CCSCalibrationTW``
        fitted TWIMS CCS calibration instance
    figname : ``str``, optional
        if not provided, do not save the figure or display it for interactive viewing. if figname is "show" then
        display the figure for interactive viewing, and if figname is a path to an image file then save the image to
        that path
    """
    # ensure that the calibration has been fitted
    if calibration.opt_params is None:
        msg = 'plot_tw_ccs_calibration: calibration has not been fit yet (calibration.opt_params is None)'
        raise RuntimeError(msg)
    # get reference and fitted values
    X, y, y_fit = calibration._X, calibration._y, calibration._y_fit
    ccs_ref = calibration.ref_ccs
    ccs_cal = calibration.calibrated_ccs(calibration.mz, calibration.arrival_time)
    percent_error = 100. * (ccs_cal - ccs_ref) / ccs_ref
    # get a text description of the fit parameters
    text_format = {
        'linear': "Y = {:.2e} * X {:+.2e}",
        'quadratic': "Y = {:.2e} * X^2 {:+.2e} * X {:+.2e}",
        'power1': "Y = {:.2e} {:+.2e} * X^{:.3f}",
        'power2': "Y = {:.2e} * (X {:+5.2e})^{:.3f}",
    }
    param_text = text_format[calibration.fit_func].format(*calibration.opt_params)
    fig, axs = plt.subplots(nrows=2, figsize=(3.33, 4), gridspec_kw={'height_ratios': [2.5, 1]})
    xlab = 'arrival time (corrected)' if calibration.correct_dt else 'arrival time'
    axs[0].plot(X, y, 'b.', ms=5)
    axs[0].plot(X, y_fit, 'b--', lw=0.75)
    axs[0].set_ylabel('CCS (corrected)' if calibration.correct_ccs else 'CCS')
    axs[0].text(0.05, 1., param_text, color='b', transform=axs[0].transAxes)
    axs[1].plot(X, percent_error, 'b.', ms=5)
    axs[1].axhline(0, ls='--', lw=0.5, color='k')
    axs[1].set_xlabel(xlab)
    axs[1].set_ylabel('CCS % error')
    for ax in axs:
        for d in ['top', 'right']:
            ax.spines[d].set_visible(False)
    if figname is not None:
        plt.tight_layout()
        if figname == 'show':
            plt.show()
        else:
            plt.savefig(figname, dpi=350, bbox_inches='tight')


