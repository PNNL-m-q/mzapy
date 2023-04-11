"""
mzapy/__init__.py

Dylan Ross (dylan.ross@pnnl.gov)
Joon-Yong Lee (junyoni@gmail.com)

    The MZA class is defined here. MZA is the main high-level interface for extracting MS data from MZA files
"""


# mza_version.major_version.minor_version
# mza_version is kept in lockstep with release of MZA format
__version__ = '1.6.3'


import queue
import threading
import time
import os
import pickle

import h5py
import hdf5plugin
import numpy as np
import pandas as pd

from mzapy._config import _MZA_VERSIONS_SUPPORTED,  _ALL_METADATA_HEADERS, _CACHE_METADATA_HEADERS
from mzapy._util import _UpdatingProgressBar


class MZA():
    """
    object for accessing MS data from HDF5 formatted files
    
    Attributes
    ----------
    mza_version : ``str``
        version of the underlying .mza format
    h5_file : ``str``
        file name (and optionally path to) the source HDF5 file
    h5 : ``h5py.File``
        h5py File object for extracting information from the HDF5 file
    mz_full : numpy.ndarray(float)
        full m/z array
    min_mz : ``float``
        minimum m/z value
    max_mz : ``float``
        maximum m/z value
    idx : ``numpy.ndarray(int)``
        indices of all scans
    dt : ``numpy.ndarray(float)``
        drift time of all scans
    min_dt : ``float``
        minimum drift time
    max_dt : ``float``
        maximum drift time
    rt : ``numpy.ndarray(float)``
        retention time of all scans
    min_rt : ``float``
        minimum retention time
    max_rt : ``float``
        maximum retention time
    imb : ``numpy.array(int)``
        ion mobility bin of all scans
    mslvl : ``numpy.ndarray(int)``
        MS level of all scans
    ms1_frame2frame_idx : ``dict(:)``
        mapping ?
    frame_idx2ms1_frame : ``dict(:)``
        mapping ?
    """

    def __init__(self, h5_file, cache_metadata_headers=_CACHE_METADATA_HEADERS, ms1lvl=1, io_threads=8, 
                 cache_scan_data=False, mza_version='new'):
        """
        init a new MZA instance using the path to the source HDF5 file

        Parameters
        ----------
        h5_file : ``str``
            file name (and optionally path to) the source HDF5 file
        cache_metadata_headers : ``list(str)`` or ``str``, default=_CACHE_METADATA_HEADERS
            specify which metadata headers to cache in memory for faster access, [] to cache none, 'all' to cache all,
            by default the set defined in _config._CACHE_METADATA_HEADERS is used
        ms1lvl : ``int``, default=1
            mslvl value corresponding to MS1 data (1 if MSMS data is present, 0 other times)
        io_threads : ``int``, default=8
            number of threads to use for performing IO tasks
        cache_scan_data : ``bool``, default=False
            whether to cache extracted scan data for faster subsequent access
        mza_version : ``str``, default='new'
            temporary measure for indicating whether the the scan indexing needs to account for partitioned
            scan data ('new') or not ('old'). Again, this is only temporary as at some point the mza version
            will be encoded as metadata into the file and this accommodation can be made automatically.
        """
        # validate and store the mza version
        if mza_version not in _MZA_VERSIONS_SUPPORTED:
            msg = 'MZA: __init__: this version of MZA does not support .mza format version: {}'
            raise ValueError(msg.format(mza_version))
        self.mza_version = mza_version
        self.h5_file = h5_file
        self.h5 = h5py.File(self.h5_file, 'r')
        # preload some of the metadata into memory for much faster access later on
        self._metadata = {}
        # check if cache_metadata_headers was set to 'all'
        cache_metadata_headers = _ALL_METADATA_HEADERS if cache_metadata_headers == 'all' else cache_metadata_headers
        for hdr in cache_metadata_headers:
            self._metadata[hdr] = self.h5['Metadata'][hdr][()]
        # deal with changes to mza format 'new' version
        if mza_version == 'new':
            # cache MzaPath metadata header if mza_version is new and cache_metadata_headers is the default or 'all'
            if cache_metadata_headers in [_CACHE_METADATA_HEADERS, 'all']:
                self._metadata['MzaPath'] = self.h5['Metadata']['MzaPath'][()]
            # map scan indices to scan partitions
            self._idx_to_path = {}
            for i, p in zip(self.metadata('Scan'), self.metadata('MzaPath')):
                self._idx_to_path[i] = p.decode()
        # preload the full m/z array into memory, it potentially gets indexed quite a few times
        self.mz_full = self.h5['Full_mz_array'][()]
        self.min_mz, self.max_mz = min(self.mz_full), max(self.mz_full)
        # set mappings between ims frame and ms1 frame
        self._ms1lvl = ms1lvl
        self._get_ms1_frames()
        # set up for threaded IO
        self._io_threads = io_threads
        self._setup_threaded_io()
        # set up scan cache mapping scan index to scan data, None means no caching is performed
        self._scan_cache = {} if cache_scan_data else None
        if self._scan_cache is not None:
            # if the cache file is not found, ignore and just start with {}
            self.load_scan_cache(ignore_no_cache_file=True)

    @property
    def _scan_cache_file(self):
        """ default filename for scan cache file, data file with extra stuff added at the end """
        return self.h5_file + '.scan_cache.pkl'

    def load_scan_cache(self, scan_cache_file=None, ignore_no_cache_file=False):
        """
        load the scan cache from file

        Parameters
        ----------
        scan_cache_file : ``str``, optional
            if provided, override the default scan cache file name
        ignore_no_cache_file : ``bool``, default=False
            do not raise an exception if the cache file is not found, just silently use {} as the scan cache
        """
        scan_cache_file = self._scan_cache_file if scan_cache_file is None else scan_cache_file
        if not os.path.isfile(scan_cache_file):
            if ignore_no_cache_file:
                self._scan_cache = {}
            else:
                msg = 'MZA.load_scan_cache: scan cache file {} not found'
                raise RuntimeError(msg.format(scan_cache_file))
        else:
            with open(scan_cache_file, 'rb') as pf:
                self._scan_cache = pickle.load(pf)

    def save_scan_cache(self, scan_cache_file=None):
        """
        saves the scan cache to file for fast loading later

        Parameters
        ----------
        scan_cache_file : ``str``, optional
            if provided, override the default scan cache file name
        """
        if self._scan_cache is None:
            msg = 'MZA.save_scan_cache: there is no scan cache to save'
            raise RuntimeError(msg)
        scan_cache_file = self._scan_cache_file if scan_cache_file is None else scan_cache_file
        with open(scan_cache_file, 'wb') as pf:
            pickle.dump(self._scan_cache, pf)

    def _io_worker_func(self):
        """
        IO worker thread function. Blocks waiting for self._io_q_in to give it a scan index and other input information,
        and when it gets one it extracts the mzbin or intensities array from the file and places those along with the 
        index as a tuple into self._io_q_out. If the thread recieves None from the input queue, breaks the main loop
        """
        while True:
            inputs = self._io_q_in.get()
            if inputs is None:
                break
            scan_idx, component, bin_idx_min, bin_idx_max = inputs
            if component == 'mzbins':
                data = self._read_scan_mzbins(scan_idx, bin_idx_min, bin_idx_max)
            elif component == 'intensities':
                data = self._read_scan_intensity(scan_idx, bin_idx_min, bin_idx_max)
            self._io_q_out.put((scan_idx, component, data))
            self._io_q_in.task_done()

    def _setup_threaded_io(self):
        """
        Sets up input/output Queues and a list of IO workers. After this method runs, the worker threads will sit and
        block waiting for scan indices to become available from self._io_q_in. Whenever a worker thread gets a scan
        index, it will extract the data for that scan and put it (along with the index) into self._io_q_out. Workers
        will stop when they recieve a None from the input Queue
        """
        # setup input/output queues, input takes scan indices, output recieves tuples of scan index, mzbins, intensities
        self._io_q_in = queue.Queue()
        self._io_q_out = queue.Queue()
        # initialize some workers and start them up
        self._io_workers = [threading.Thread(target=self._io_worker_func) for _ in range(self._io_threads)]
        for wkr in self._io_workers:
            wkr.start()

    def _stop_io_workers(self):
        """
        stops all IO worker threads (which otherwise would be blocking waiting for the input queue to give them some 
        input data) by pushing a bunch of Nones into the input queue. This is a shutdown task.
        """
        for _ in range(self._io_threads):
            self._io_q_in.put(None)

    def close(self):
        """
        closes the underlying h5py File object, trying to read data after closing will cause errors
        """
        self.h5.close()
        # TODO (Dylan Ross): What to do if the input queue has unprocessed inputs in it? pushing Nones will do nothing
        #                    until the threads finish what they are working on and make it to the Nones in the input
        #                    queue. In order to close gracefully, even when processing has been interrupted for some 
        #                    reason, this situation needs to be detected and dealt with here. 
        # stop io worker threads
        self._stop_io_workers()
        # if scan cache is turned on, save it to file on close
        if self._scan_cache is not None:
            self.save_scan_cache()

    def metadata(self, header):
        """
        retrieve the array of metadata from the specified header

        Parameters
        ----------
        header : ``str``
            specify which metadata header to fetch

        Returns
        -------
        metadata_column : ``numpy.ndarray(?)``
            array of metadata (could be type ``float``, ``int``, or ``bytes`` depending on the header)
        """
        if header in self._metadata:
            return self._metadata[header]
        else:
            return self.h5['Metadata'][header][()]

    @property
    def idx(self):
        """ indices of all scans """
        return self.metadata('Scan')

    @property
    def dt(self):
        """ drift time of all scans """
        return self.metadata('IonMobilityTime')

    @property
    def min_dt(self):
        """ minimum drift time"""
        return min(self.metadata('IonMobilityTime'))

    @property
    def max_dt(self):
        """ maximum drift time"""
        return max(self.metadata('IonMobilityTime'))
    
    @property
    def rt(self):
        """ retention time of all scans """
        return self.metadata('RetentionTime')

    @property
    def min_rt(self):
        """ minimum retention time"""
        return min(self.metadata('RetentionTime'))

    @property
    def max_rt(self):
        """ maximum retention time"""
        return max(self.metadata('RetentionTime'))

    @property
    def imf(self):
        """ ion mobility frame of all scans """
        return self.metadata('IonMobilityFrame')

    @property
    def imb(self):
        """ ion mobility bin of all scans """
        return self.metadata('IonMobilityBin')

    @property
    def mslvl(self):
        """ MS level of all scans """
        return self.metadata('MSLevel')

    def _closest_mzbin(self, mz, direction):
        """
        finds the closes mzbin value to an input m/z.
        when used for computing bounds it is useful to find the closest mzbin that is above or below a value, which 
        is controlled by the direction kwarg. By default, direction is 'any' meaning the closes mzbin is returned
        regardless of direction. 

        Parameters
        ----------
        mz : ``float``
            target m/z
        direction : ``str``, default='any'
            specify whether the closest mzbin above or below the target value should be returned, 'any' to ignore 
            direction and just return the closest mz bin
        """
        idx_any = dmz = np.abs(mz - self.mz_full).argmin()
        if direction == 'any':
            return idx_any
        if direction == 'above':
            if mz > self.max_mz:
                'MZA._closest_mzbin: target m/z {} is greater than maximum m/z value {}'
                raise ValueError(msg.format(mz, self.max_mz))
            return idx_any if self.mz_full[idx_any] > mz else idx_any + 1
        elif direction == 'below':
            if mz < self.min_mz:
                'MZA._closest_mzbin: target m/z {} is less than minimum m/z value {}'
                raise ValueError(msg.format(mz, self.min_mz))
            return idx_any if self.mz_full[idx_any] < mz else idx_any - 1
        msg = 'MZA._closest_mzbin: direction kwarg must be "above", "below" or "any" (was: "{}")'
        raise ValueError(msg.format(direction))
    
    def _read_scan_intensity(self, scan_idx, bin_idx_min, bin_idx_max):
        """
        reads the intensities for a specified scan index, min/max bin indices can be provided to specify that only part
        of the scan data should be extracted. If both are set to None, then the full scan data is extracted

        Parameters
        ----------
        scan_idx : ``int``
            scan index to read
        bin_idx_min : ``int`` or ``None``
            index of minimum bin to extract from scan data
        bin_idx_max : ``int`` or ``None``
            index of maximum bin to extract from scan data

        Returns
        -------
        intensity : ``numpy.ndarray(int)``
            array of intensities from specified scan
        """
        path = '{}'.format(self._idx_to_path[scan_idx]) if self.mza_version == 'new' else ''
        full_index = 'Arrays_intensity{}/{}'.format(path, scan_idx)
        if bin_idx_min is None and bin_idx_max is None:
            # full scan data
            return self.h5[full_index][()]
        else:
            # indexed scan data
            return self.h5[full_index][bin_idx_min:bin_idx_max]
    
    def _read_scan_mzbins(self, scan_idx, bin_idx_min, bin_idx_max):
        """
        reads the mzbins for a specified scan index, min/max bin indices can be provided to specify that only part of
        the scan data should be extracted. If both are set to None, then the full scan data is extracted

        Parameters
        ----------
        scan_idx : ``int``
            scan index to read
        bin_idx_min : ``int`` or ``None``
            index of minimum bin to extract from scan data
        bin_idx_max : ``int`` or ``None``
            index of maximum bin to extract from scan data

        Returns
        -------
        mzbin : numpy.ndarray(int)
            array of mzbins from specified scan
        """
        path = '{}'.format(self._idx_to_path[scan_idx]) if self.mza_version == 'new' else ''
        full_index = 'Arrays_mzbin{}/{}'.format(path, scan_idx)
        if bin_idx_min is None and bin_idx_max is None:
            # full scan data
            return self.h5[full_index][()]
        else:
            # indexed scan data
            return self.h5[full_index][bin_idx_min:bin_idx_max]

    def _get_ms1_frames(self):
        """
        initializes class attributes that map between MS1 frames and IonMobilityFrame:
            self.ms1_frame2frame_idx
            self.frame_idx2ms1_frame 

        TODO (Dylan Ross) I don't really understand what this is doing yet...
                          It seems like it just ends up mapping N : N+1 or N : N-1 ?
                          Once I figure out what this is actually achieving I will clean it up
        """
        meta = self.h5["Metadata"]["Scan", "MSLevel", 'IonMobilityFrame', 'RetentionTime']
        meta = meta[meta["MSLevel"] == self._ms1lvl]
        tdf = pd.DataFrame(meta[:], columns=["Scan", "MSLevel", 'IonMobilityFrame', 'RetentionTime'])
        frame_rt_list = tdf[['IonMobilityFrame', 'RetentionTime']].drop_duplicates().to_dict('records')
        self.ms1_frame2frame_idx = {}
        self.frame_idx2ms1_frame = {}
        for idx, row in enumerate(frame_rt_list):
            self.ms1_frame2frame_idx[idx] = row['IonMobilityFrame']
            self.frame_idx2ms1_frame[row['IonMobilityFrame']] = idx

    def _read_multi_scan_data(self, scan_idxs, verbose):
        """
        uses multithreaded IO to read scan data from the designated scan indices

        Parameters
        ----------
        scan_idxs : ``list(int)``
            scan indices for scan data
        verbose : ``bool``
            whether to print progress information

        Returns
        -------
        scan_data : ``dict(int:list(numpy.ndarray(int)))``
            dictionary mapping scan index to [mzbins, intensities] as numpy.ndarrays
        """
        scan_data = {}
        n_cached = 0
        n_uncached = 0
        for idx in scan_idxs:
            if self._scan_cache is not None and idx in self._scan_cache:
                scan_data[idx] = self._scan_cache[idx]
                n_cached += 1
            else:
                # scan indices into the input IO queue, IO threads will immediately start processing
                self._io_q_in.put((idx, 'mzbins', None, None))
                self._io_q_in.put((idx, 'intensities', None, None))
                n_uncached += 1
        # track progress
        if verbose:
            total_scans = n_uncached
            pbar = _UpdatingProgressBar(levels=40, track_time_elapsed=True)
            msg = 'found {} cached scans, extracting {} scans ...'
            print(msg.format(n_cached, total_scans), flush=True)
            progress = 100. * (1. - (self._io_q_in.qsize() / (total_scans * 2))) if n_uncached > 0 else 99 
            pbar.update(progress)
            while progress < 95:
                time.sleep(3)  # 3 seconds between trying to update the progress bar
                progress = 100. * (1. - (self._io_q_in.qsize() / (total_scans * 2)))
                pbar.update(progress)
        # wait for all threads to finish extracting data
        self._io_q_in.join()
        if verbose and not pbar.finished: 
            pbar.complete()
        # retrieve scan data from output IO queue, coordinate separate mzbin/intensities
        while not self._io_q_out.empty():
            idx, component, data = self._io_q_out.get()
            if idx not in scan_data:
                scan_data[idx] = [None, None]
            if component == 'mzbins':
                scan_data[idx][0] = data
            if component == 'intensities':
                scan_data[idx][1] = data
        # add any new scan data to cache if scan caching is turned on
        if self._scan_cache is not None:
            for idx in scan_data:
                if idx not in self._scan_cache:
                    self._scan_cache[idx] = scan_data[idx]
        return scan_data

    def collect_ms1_df_by_rt(self, rt_min, rt_max, mz_bounds=None, verbose=False):
        """
        collects MS1 data within RT ranges, IM dimension is collapsed/ignored

        Parameters
        ----------
        rt_min : ``float``
            lower bound of RT window to select data from
        rt_max : ``float``
            upper bound of RT window to select data from
        mz_bounds : ``tuple(float, float)``, optional
            mz_min, mz_max
        verbose : ``bool``, default=False
            print information about the progress

        Returns
        -------
        data : ``pandas.DataFrame``
            data frame with columns mzbin, mz, intensity, rt, frame
        """
        # determine the indices to select
        sel = (self.mslvl == self._ms1lvl) & (self.imb == 0) & (rt_min <= self.rt) & (self.rt <= rt_max)
        multi_scan_data = self._read_multi_scan_data(self.idx[sel], verbose) 
        if mz_bounds is not None:
            mz_min, mz_max = mz_bounds
        data = []
        for idx, msl, rt, frame in zip(self.idx[sel], self.mslvl[sel], self.rt[sel], 
                                       self.metadata('IonMobilityFrame')[sel]):
            ms1_frame = self.frame_idx2ms1_frame[frame]
            mzbins, intensities = multi_scan_data[idx]
            for mzbin, intensity in zip(mzbins, intensities):
                if mz_bounds is None:
                    data.append([mzbin, self.mz_full[mzbin], intensity, rt, ms1_frame])
                else:
                    mz = self.mz_full[mzbin]
                    if mz >= mz_min and mz <= mz_max:
                        data.append([mzbin, mz, intensity, rt, ms1_frame])
        return pd.DataFrame(data, columns=['mzbin', 'mz', 'intensity', 'rt', 'frame'])

    def collect_ms1_df_by_rt_dt(self, rt_min, rt_max, dt_min, dt_max, mz_bounds=None, verbose=False):
        """
        collects MS1 data within RT and DT ranges

        Parameters
        ----------
        rt_max : ``float``
            lower bound of RT window to select data from
        rt_max : ``float``
            upper bound of RT window to select data from
        dt_min : ``float``
            lower bound of DT window to select data from
        dt_max : ``float``
            upper bound of DT window to select data from
        mz_bounds : ``tuple(float, float)``, optional
            mz_min, mz_max
        verbose : ``bool``, default=False
            print information about the progress

        Returns
        -------
        data : ``pandas.DataFrame``
            data frame with columns mzbin, mz, intensity, rt, dt, frame
        """
        sel = (self.mslvl == self._ms1lvl) & (self.imb != 0) & (rt_min <= self.rt) & \
              (self.rt <= rt_max) & (dt_min <= self.dt) & (self.dt <= dt_max)
        multi_scan_data = self._read_multi_scan_data(self.idx[sel], verbose)
        if mz_bounds is not None:
            mz_min, mz_max = mz_bounds
        data = []  
        for idx, msl, rt, dt, frame in zip(self.idx[sel], self.mslvl[sel], self.rt[sel], self.dt[sel],
                                           self.metadata('IonMobilityFrame')[sel]):
            mzbins, intensities = multi_scan_data[idx]
            for mzbin, intensity in zip(mzbins, intensities):
                if mz_bounds is None:
                    data.append([mzbin, self.mz_full[mzbin], intensity, rt, dt, self.frame_idx2ms1_frame[frame]])
                else:
                    mz = self.mz_full[mzbin]
                    if mz >= mz_min and mz <= mz_max:
                        data.append([mzbin, self.mz_full[mzbin], intensity, rt, dt, self.frame_idx2ms1_frame[frame]])
        return pd.DataFrame(data, columns=['mzbin', 'mz', 'intensity', 'rt', 'dt', 'frame'])

    def collect_ms1_arrays_by_rt(self, rt_min, rt_max, mz_bounds=None):
        """
        loads MS1 spectrum (m/z, intensity) as arrays for m/z and RT range, *ignoring DT if present*

        Parameters
        ----------
        rt_min : ``float``
            lower RT bound
        rt_max : ``float``
            upper RT bound
        mz_bounds : ``tuple(float)``, optional
            (lower, upper) m/z bounds, filters data after extraction so no effect on extraction time

        Returns
        -------
        ms1_mz : ``numpy.ndarray(float)``
            m/z component of MS1 spectrum
        ms1_int : ``numpy.ndarray(int)``
            intensity component of MS1 spectrum
        """
        ms1_df = self.collect_ms1_df_by_rt(rt_min, rt_max, mz_bounds=mz_bounds)
        df = ms1_df.groupby('mz').intensity.sum()
        ms1_mz, ms1_int = df.index.to_numpy(), df.to_numpy()
        return ms1_mz, ms1_int

    def collect_ms1_arrays_by_dt(self, dt_min, dt_max, mz_bounds=None):
        """
        loads MS1 spectrum (m/z, intensity) as arrays for m/z and DT range

        Parameters
        ----------
        dt_min : ``float``
            lower DT bound
        dt_max : ``float``
            upper DT bound
        mz_bounds : ``tuple(float)``, optional
            (lower, upper) m/z bounds, filters data after extraction so no effect on extraction time

        Returns
        -------
        ms1_mz : ``numpy.ndarray(float)``
            m/z component of MS1 spectrum
        ms1_int : ``numpy.ndarray(int)``
            intensity component of MS1 spectrum
        """
        ms1_df = self.collect_ms1_df_by_rt_dt(self.min_rt, self.max_rt, dt_min, dt_max, mz_bounds=mz_bounds)
        df = ms1_df.groupby('mz').intensity.sum()
        ms1_mz, ms1_int = df.index.to_numpy(), df.to_numpy()
        return ms1_mz, ms1_int

    def collect_ms2_df_by_rt(self, rt_min, rt_max, mz_bounds=None, verbose=False):
        """
        collects MS2 data within RT ranges, IM dimension is collapsed/ignored

        Parameters
        ----------
        rt_min : ``float``
            lower bound of RT window to select data from
        rt_max : ``float``
            upper bound of RT window to select data from
        mz_bounds : ``tuple(float, float)``, optional
            mz_min, mz_max
        verbose : ``bool``, default=False
            print information about the progress

        Returns
        -------
        data : ``pandas.DataFrame``
            data frame with columns mzbin, mz, intensity, rt
        """
        # determine the indices to select
        sel = (self.mslvl == 2) & (self.imb == 0) & (rt_min <= self.rt) & (self.rt <= rt_max)
        multi_scan_data = self._read_multi_scan_data(self.idx[sel], verbose) 
        if mz_bounds is not None:
            mz_min, mz_max = mz_bounds
        data = []
        for idx, msl, rt, frame in zip(self.idx[sel], self.mslvl[sel], self.rt[sel], 
                                       self.metadata('IonMobilityFrame')[sel]):
            mzbins, intensities = multi_scan_data[idx]
            for mzbin, intensity in zip(mzbins, intensities):
                if mz_bounds is None:
                    data.append([mzbin, self.mz_full[mzbin], intensity, rt])
                else:
                    mz = self.mz_full[mzbin]
                    if mz >= mz_min and mz <= mz_max:
                        data.append([mzbin, mz, intensity, rt])
        return pd.DataFrame(data, columns=['mzbin', 'mz', 'intensity', 'rt'])

    def collect_ms2_df_by_rt_dt(self, rt_min, rt_max, dt_min, dt_max, mz_bounds=None, verbose=False):
        """
        collects MS2 data within RT and DT ranges

        Parameters
        ----------
        rt_max : ``float``
            lower bound of RT window to select data from
        rt_max : ``float``
            upper bound of RT window to select data from
        dt_min : ``float``
            lower bound of DT window to select data from
        dt_max : ``float``
            upper bound of DT window to select data from
        mz_bounds : ``tuple(float, float)``, optional
            mz_min, mz_max
        verbose : ``bool``, default=False
            print information about the progress

        Returns
        -------
        data : ``pandas.DataFrame``
            data frame with columns mzbin, mz, intensity, rt, dt
        """
        sel = (self.mslvl == 2) & (self.imb != 0) & (rt_min <= self.rt) & \
              (self.rt <= rt_max) & (dt_min <= self.dt) & (self.dt <= dt_max)
        multi_scan_data = self._read_multi_scan_data(self.idx[sel], verbose)
        if mz_bounds is not None:
            mz_min, mz_max = mz_bounds
        data = []  
        for idx, msl, rt, dt, frame in zip(self.idx[sel], self.mslvl[sel], self.rt[sel], self.dt[sel],
                                           self.metadata('IonMobilityFrame')[sel]):
            mzbins, intensities = multi_scan_data[idx]
            for mzbin, intensity in zip(mzbins, intensities):
                if mz_bounds is None:
                    data.append([mzbin, self.mz_full[mzbin], intensity, rt, dt])
                else:
                    mz = self.mz_full[mzbin]
                    if mz >= mz_min and mz <= mz_max:
                        data.append([mzbin, self.mz_full[mzbin], intensity, rt, dt])
        return pd.DataFrame(data, columns=['mzbin', 'mz', 'intensity', 'rt', 'dt'])

    def collect_ms2_arrays_by_rt(self, rt_min, rt_max, mz_bounds=None):
        """
        loads MS2 spectrum (m/z, intensity) as arrays for m/z and RT range, *ignoring DT if present*

        Parameters
        ----------
        rt_min : ``float``
            lower RT bound
        rt_max : ``float``
            upper RT bound
        mz_bounds : ``tuple(float)``, optional
            (lower, upper) m/z bounds, filters data after extraction so no effect on extraction time

        Returns
        -------
        ms1_mz : ``numpy.ndarray(float)``
            m/z component of MS1 spectrum
        ms1_int : ``numpy.ndarray(int)``
            intensity component of MS1 spectrum
        """
        ms1_df = self.collect_ms2_df_by_rt(rt_min, rt_max, mz_bounds=mz_bounds)
        df = ms1_df.groupby('mz').intensity.sum()
        ms1_mz, ms1_int = df.index.to_numpy(), df.to_numpy()
        return ms1_mz, ms1_int

    def collect_ms2_arrays_by_rt_dt(self, rt_min, rt_max, dt_min, dt_max, mz_bounds=None, verbose=False):
        """
        loads MS2 spectrum (m/z, intensity) as arrays for m/z, RT, and DT range

        Parameters
        ----------
        rt_min : ``float``
            lower RT bound
        rt_max : ``float``
            upper RT bound
        dt_min : ``float``
            lower bound of DT window to select data from
        dt_max : ``float``
            upper bound of DT window to select data from
        mz_bounds : ``tuple(float)``, optional
            (lower, upper) m/z bounds, filters data after extraction so no effect on extraction time
        verbose : ``bool``, default=False
            print information about the progress

        Returns
        -------
        ms1_mz : ``numpy.ndarray(float)``
            m/z component of MS1 spectrum
        ms1_int : ``numpy.ndarray(int)``
            intensity component of MS1 spectrum
        """
        ms1_df = self.collect_ms2_df_by_rt_dt(rt_min, rt_max, dt_min, dt_max, 
                                              mz_bounds=mz_bounds, verbose=verbose)
        df = ms1_df.groupby('mz').intensity.sum()
        ms1_mz, ms1_int = df.index.to_numpy(), df.to_numpy()
        return ms1_mz, ms1_int

    def collect_xic_arrays_by_mz(self, mz_min, mz_max, rt_bounds=None, mslvl=1, verbose=False):
        """
        loads XIC (retention time, intensity) as arrays for m/z range, *ignoring DT if present*

        Parameters
        ----------
        mz_min : ``float``
            lower m/z bound
        mz_max : ``float``
            upper m/z bound
        rt_bounds : ``tuple(float, float)``, optional
            (lower, upper) RT bounds
        mslvl : ``int``, default=1
            MS level to select from
        verbose : ``bool``, default=False
            print information about the progress

        Returns
        -------
        xic_rt : ``numpy.ndarray(float)``
            retention time component of XIC
        xic_int : ``numpy.ndarray(int)``
            intensity component of XIC
        """
        if rt_bounds is not None:
            rt_min, rt_max = rt_bounds
        else:
            rt_min, rt_max = self.min_rt, self.max_rt
        if mslvl == 1: 
            df = self.collect_ms1_df_by_rt(rt_min, rt_max, mz_bounds=(mz_min, mz_max), verbose=verbose).groupby('rt').intensity.sum()
        elif mslvl ==2 :
            df = self.collect_ms2_df_by_rt(rt_min, rt_max, mz_bounds=(mz_min, mz_max), verbose=verbose).groupby('rt').intensity.sum()
        xic_rt, xic_int = df.index.to_numpy(), df.to_numpy()
        return xic_rt, xic_int

    def collect_ms1_arrays_by_rt_dt(self, rt_min, rt_max, dt_min, dt_max, mz_bounds=None, verbose=False):
        """
        loads MS1 spectrum (m/z, intensity) as arrays for m/z, RT, and DT range

        Parameters
        ----------
        rt_min : ``float``
            lower RT bound
        rt_max : ``float``
            upper RT bound
        dt_min : ``float``
            lower bound of DT window to select data from
        dt_max : ``float``
            upper bound of DT window to select data from
        mz_bounds : ``tuple(float)``, optional
            (lower, upper) m/z bounds, filters data after extraction so no effect on extraction time
        verbose : ``bool``, default=False
            print information about the progress

        Returns
        -------
        ms1_mz : ``numpy.ndarray(float)``
            m/z component of MS1 spectrum
        ms1_int : ``numpy.ndarray(int)``
            intensity component of MS1 spectrum
        """
        ms1_df = self.collect_ms1_df_by_rt_dt(rt_min, rt_max, dt_min, dt_max, 
                                              mz_bounds=mz_bounds, verbose=verbose)
        df = ms1_df.groupby('mz').intensity.sum()
        ms1_mz, ms1_int = df.index.to_numpy(), df.to_numpy()
        return ms1_mz, ms1_int

    def collect_atd_arrays_by_rt_mz(self, mz_min, mz_max, rt_min, rt_max, dt_bounds=None, mslvl=1, verbose=False):
        """
        loads ATD (dt, intensity) as arrays for target mass within an RT window

        Parameters
        ----------
        mz_min : ``float``
            lower m/z bound
        mz_max : ``float``
            upper m/z bound
        rt_min : ``float``
            lower RT bound
        rt_max : ``float``
            upper RT bound
        dt_bounds : ``tuple(float)``, optional
            (lower, upper) drift bounds, tightening DT bounds around area of interest reduces extraction time
        mslvl : ``int``, default=1
            MS level to select from 
        verbose : ``bool``, default=False
            print information about the progress

        Returns
        -------
        atd_dt : ``numpy.ndarray(float)``
            drift time component of ATD
        atd_int : ``numpy.ndarray(int)``
            intensity component of ATD
        """
        if dt_bounds is not None:
            dt_min, dt_max = dt_bounds
        else:
            dt_min, dt_max = self.min_dt, self.max_dt
        if mslvl == 1:
            df = self.collect_ms1_df_by_rt_dt(rt_min, rt_max, dt_min, dt_max, 
                                            mz_bounds=(mz_min, mz_max), verbose=verbose).groupby('dt').intensity.sum()
        elif mslvl == 2:
            df = self.collect_ms2_df_by_rt_dt(rt_min, rt_max, dt_min, dt_max, 
                                            mz_bounds=(mz_min, mz_max), verbose=verbose).groupby('dt').intensity.sum()
        atd_dt, atd_int = df.index.to_numpy(), df.to_numpy()
        return atd_dt, atd_int

    def collect_xic_arrays_by_mz_dt(self, mz_min, mz_max, dt_min, dt_max, rt_bounds=None, mslvl=1, verbose=False):
        """
        loads XIC (RT, intensity) as arrays for target mass within a DT window

        Parameters
        ----------
        mz_min : ``float``
            lower m/z bound
        mz_max : ``float``
            upper m/z bound
        dt_min : ``float``
            lower DT bound
        dt_max : ``float``
            upper DT bound
        rt_bounds : ``tuple(float, float)``, optional
            (lower, upper) RT bounds
        mslvl : ``int``, default=1
            MS level to select from
        verbose : ``bool``, default=False
            print information about the progress

        Returns
        -------
        xic_rt : ``numpy.ndarray(float)``
            retention time component of XIC
        xic_int : ``numpy.ndarray(int)``
            intensity component of XIC
        """
        if rt_bounds is not None:
            rt_min, rt_max = rt_bounds
        else:
            rt_min, rt_max = self.min_rt, self.max_rt
        if mslvl == 1:
            df = self.collect_ms1_df_by_rt_dt(rt_min, rt_max, dt_min, dt_max, 
                                              mz_bounds=(mz_min, mz_max), verbose=verbose).groupby('rt').intensity.sum()
        elif mslvl == 2:
            df = self.collect_ms2_df_by_rt_dt(rt_min, rt_max, dt_min, dt_max, 
                                              mz_bounds=(mz_min, mz_max), verbose=verbose).groupby('rt').intensity.sum()
        xic_rt, xic_int = df.index.to_numpy(), df.to_numpy()
        return xic_rt, xic_int

    def collect_dtmz_arrays_by_rt(self, rt_min, rt_max, dt_bounds=None, mz_bounds=None, verbose=False):
        """ 
        loads DTMZ data (dt, mz, intensity) as arrays within a target RT range, using optional m/z and DT bounds 
        
        Parameters
        ----------
        rt_min : ``float``
            lower RT bound
        rt_max : ``float``
            upper RT bound
        dt_bounds : ``tuple(float)``, optional
            (lower, upper) drift bounds, tightening DT bounds around area of interest reduces extraction time
        mz_bounds : ``tuple(float)``, optional
            (lower, upper) m/z bounds, filters data after extraction so no effect on extraction time
        verbose : ``bool``, default=False
            print information about the progress

        Returns
        -------
        dtmz_dt : ``numpy.ndarray(float)``
            drift time component of DTMZ data
        dtmz_mz : ``numpy.ndarray(float)``
            m/z component of DTMZ data
        dtmz_int : ``numpy.ndarray(int)``
            intensity component of DTMZ data
        """
        if dt_bounds is not None:
            dt_min, dt_max = dt_bounds
        else:
            dt_min, dt_max = self.min_dt, self.max_dt
        ms1df = self.collect_ms1_df_by_rt_dt(rt_min, rt_max, dt_min, dt_max, mz_bounds=mz_bounds, verbose=verbose)
        df2 = ms1df.groupby(['dt', 'mz']).intensity.sum()
        dtmz_dt = df2.index.get_level_values('dt').to_numpy()
        dtmz_mz = df2.index.get_level_values('mz').to_numpy()
        dtmz_i = df2.to_numpy()
        return dtmz_dt, dtmz_mz, dtmz_i

    def collect_rtmz_arrays_by_dt(self, dt_min, dt_max, rt_bounds=None, mz_bounds=None, verbose=False):
        """ 
        loads RTMZ data (rt, mz, intensity) as arrays within a target DT range, using optional m/z and RT bounds 
        
        Parameters
        ----------
        dt_min : ``float``
            lower DT bound
        dt_max : ``float``
            upper DT bound
        rt_bounds : ``tuple(float)``, optional
            (lower, upper) retention time bounds, tightening RT bounds around area of interest reduces extraction time
        mz_bounds : ``tuple(float)``, optional
            (lower, upper) m/z bounds, filters data after extraction so no effect on extraction time
        verbose : ``bool``, default=False
            print information about the progress

        Returns
        -------
        rtmz_rt : ``numpy.ndarray(float)``
            retention time component of RTMZ data
        rtmz_mz : ``numpy.ndarray(float)``
            m/z component of RTMZ data
        rtmz_int : ``numpy.ndarray(int)``
            intensity component of RTMZ data
        """
        if rt_bounds is not None:
            rt_min, rt_max = rt_bounds
        else:
            rt_min, rt_max = self.min_rt, self.max_rt
        ms1df = self.collect_ms1_df_by_rt_dt(rt_min, rt_max, dt_min, dt_max, mz_bounds=mz_bounds, verbose=verbose)
        df2 = ms1df.groupby(['rt', 'mz']).intensity.sum()
        rtmz_rt = df2.index.get_level_values('rt').to_numpy()
        rtmz_mz = df2.index.get_level_values('mz').to_numpy()
        rtmz_i = df2.to_numpy()
        return rtmz_rt, rtmz_mz, rtmz_i

    def collect_rtmz_arrays(self, rt_bounds=None, mz_bounds=None, verbose=False):
        """ 
        loads RTMZ data (rt, mz, intensity) as arrays ignoring/collapsing DT, using optional m/z and RT bounds 
        
        Parameters
        ----------
        rt_bounds : ``tuple(float)``, optional
            (lower, upper) retention time bounds, tightening RT bounds around area of interest reduces extraction time
        mz_bounds : ``tuple(float)``, optional
            (lower, upper) m/z bounds, filters data after extraction so no effect on extraction time
        verbose : ``bool``, default=False
            print information about the progress

        Returns
        -------
        rtmz_rt : ``numpy.ndarray(float)``
            retention time component of RTMZ data
        rtmz_mz : ``numpy.ndarray(float)``
            m/z component of RTMZ data
        rtmz_int : ``numpy.ndarray(int)``
            intensity component of RTMZ data
        """
        if rt_bounds is not None:
            rt_min, rt_max = rt_bounds
        else:
            rt_min, rt_max = self.min_rt, self.max_rt
        ms1df = self.collect_ms1_df_by_rt(rt_min, rt_max, mz_bounds=mz_bounds, verbose=verbose)
        df2 = ms1df.groupby(['rt', 'mz']).intensity.sum()
        rtmz_rt = df2.index.get_level_values('rt').to_numpy()
        rtmz_mz = df2.index.get_level_values('mz').to_numpy()
        rtmz_i = df2.to_numpy()
        return rtmz_rt, rtmz_mz, rtmz_i

    def collect_rtdt_arrays_by_mz(self, mz_min, mz_max, rt_bounds=None, dt_bounds=None, verbose=False):
        """ 
        loads RTDT data (rt, dt, intensity) as arrays within a target m/z range, using optional RT and DT bounds 
        
        Parameters
        ----------
        mz_min : ``float``
            lower m/z bound
        mz_max : ``float``
            upper m/z bound
        rt_bounds : ``tuple(float)``, optional
            (lower, upper) retention time bounds, tightening RT bounds around area of interest reduces extraction time
        dt_bounds : ``tuple(float)``, optional
            (lower, upper) drift time bounds, tightening DT bounds around area of interest reduces extraction time
        verbose : ``bool``, default=False
            print information about the progress

        Returns
        -------
        rtdt_rt : ``numpy.ndarray(float)``
            retention time component of RTDT data
        rtdt_dt : ``numpy.ndarray(float)``
            drift time component of RTDT data
        rtdt_int : ``numpy.ndarray(int)``
            intensity component of RTDT data
        """
        if rt_bounds is not None:
            rt_min, rt_max = rt_bounds
        else:
            rt_min, rt_max = self.min_rt, self.max_rt
        if dt_bounds is not None:
            dt_min, dt_max = dt_bounds
        else:
            dt_min, dt_max = self.min_dt, self.max_dt
        mz_bounds = (mz_min, mz_max)
        ms1df = self.collect_ms1_df_by_rt_dt(rt_min, rt_max, dt_min, dt_max, mz_bounds=mz_bounds, verbose=verbose)
        df2 = ms1df.groupby(['rt', 'dt']).intensity.sum()
        rtdt_rt = df2.index.get_level_values('rt').to_numpy()
        rtdt_dt = df2.index.get_level_values('dt').to_numpy()
        rtdt_i = df2.to_numpy()
        return rtdt_rt, rtdt_dt, rtdt_i

