"""
mzapy/_util.py

Dylan Ross (dylan.ross@pnnl.gov)

    internal module with general utility functions
"""


import time

import numpy as np


def _ppm_error(reference, other):
    """
    compute the error in ppm for a value relative to some reference value
    
    Paramters
    ---------
    reference : float
        reference value
    other : float
        other value

    Returns
    -------
    ppm : float
        ppm error of other relative to reference
    """
    return 1e6 * (other - reference) / reference


def _abs_ppm_error(reference, other):
    """
    compute the absolute error in ppm for a value relative to some reference value
    
    Paramters
    ---------
    reference : float
        reference value
    other : float
        other value

    Returns
    -------
    ppm : float
        ppm error of other relative to reference
    """
    return 1e6 * np.abs(other - reference) / reference


class _UpdatingProgressBar():
    """ a simple class for reporting progress """

    def __init__(self, levels=20, track_time_elapsed=False):
        """
        inits a new _UpdatingProgressBar instance with a specified number of levels

        Parameters
        ----------
        levels : int, default=20
            levels in the complete bar
        track_time_elapsed : bool, default=False
            if True, keep track of the time that elapses between calls to __init__() and complete(), report
            the elapsed time (in seconds) in the completion message
        """
        self.levels = levels
        self.progress = 0.
        self.finished = False
        self.track_time_elapsed = track_time_elapsed
        if self.track_time_elapsed:
            self.t0 = time.time()

    def _print_bar(self):
        """
        prints a progress bar with the current value of self.progress. Uses carriage return to overwrite previous
        progress bar 
        """
        progress_int = int(self.levels * self.progress / 100.) 
        s = '\rprogress [' + ''.join(['=' for _ in range(progress_int)]) + ''.join([' ' for _ in range(self.levels - progress_int)]) + ']'
        s += ' {:6.2f} %'.format(self.progress)
        print(s, end=' ', flush=True)

    def update(self, percent):
        """
        updates the state of the progress bar with a specified percentage completion. 
        if percent goes above 100, call self.complete to print the final bar

        automatically re-prints the bar after update
        
        Parameters
        ----------
        percent : float
            new progress level
        """
        self.progress = percent
        if int(self.progress) < 0:
            # no negative progress allowed
            self.progress = 0
        if int(self.progress) >= 100: 
            self.complete()
        else:
            self._print_bar()

    def increment(self, percent):
        """
        increments the state of the progress bar by a specified percentage completion.
        
        Parameters
        ----------
        percent : float or int
            new progress level
        """
        self.update(self.progress + percent)

    def complete(self):
        """
        prints a full progress bar with a newline and the word "done"
        """
        self.finished = True
        self.progress = 100.
        self._print_bar()
        complete_msg = 'done'
        if self.track_time_elapsed:
            complete_msg += ' ({:.1f} s)'.format(time.time() - self.t0)
        print(complete_msg, flush=True)

