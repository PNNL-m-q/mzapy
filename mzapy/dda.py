"""
mzapy/dda.py

Dylan Ross (dylan.ross@pnnl.gov)

    Define a couple of reader objects for dealing with MZA files with DDA MS/MS data 
"""


import h5py
import hdf5plugin
import numpy as np
import pandas as pd


class MsmsReaderDda():
    """
    object for accessing DDA MSMS data from MZA (HDF5) formatted files

    The MZA files created from qTOF DDA data have some peculiarities so a
    small, purpose-built reader is better than shoehorning the main MZA object for
    this purpose

    Attributes
    ----------
    h5 : ``h5py.File``
        interface to MZA file (HDF5)
    f : ``str``
        path to MZA file
    metadata : ``pandas.DataFrame``
        metadata dataframe
    arrays_mz : ``pandas.DataFrame``
        m/z dataframe
    arrays_i : ``pandas.DataFrame``
        intensity dataframe
    ms1_scans : ``numpy.ndarray(int)``
        array of MS1 scans
    min_rt, max_rt : ``float``
        minimum/maximum retention time
    """

    def __init__(self, mza_file, drop_scans=None):
        """
        initialize the reader
    
        Parameters
        ----------
        mza_file : ``str``
            path to MZA file
        drop_scans : ``list(int)``, optional
            list of scans to drop from the file, can be None if there are not any to drop
        """
        self.h5 = h5py.File(mza_file, 'r')
        self.f = mza_file
        self.metadata = pd.DataFrame(self.h5['Metadata'][:]).set_index('Scan')
        self.arrays_mz = pd.DataFrame(self.h5['Arrays_mzbin'].items(), columns=['Scan', 'Data']).set_index('Scan')
        self.arrays_mz.index = self.arrays_mz.index.astype('int64')
        self.arrays_i = pd.DataFrame(self.h5['Arrays_intensity'].items(), columns=['Scan', 'Data']).set_index('Scan')
        self.arrays_i.index = self.arrays_i.index.astype('int64')
        if drop_scans is not None:
            #print('dropping scans:', drop_scans)
            self.metadata.drop(drop_scans, inplace=True)
            self.arrays_mz.drop(drop_scans, inplace=True)
            self.arrays_i.drop(drop_scans, inplace=True)
        self.ms1_scans = self.metadata[self.metadata['MSLevel'] == 1].index.to_numpy()
        self.ms2_scans = self.metadata[self.metadata['MSLevel'] == 2].index.to_numpy()
        rt = self.metadata[self.metadata['MSLevel'] == 1].loc[:, 'RetentionTime'].to_numpy()
        self.min_rt, self.max_rt = min(rt), max(rt)
        self.mz_full = self.h5['Full_mz_array'][()]
        self.min_mz, self.max_mz = min(self.mz_full), max(self.mz_full)
    
    def close(self):
        """
        close connection to the MZA file
        """
        self.h5.close()
        
    def get_chrom(self, mz, mz_tol, rt_bounds=None):
        """
        Select a chromatogram (MS1 only) for a target m/z with specified tolerance

        Parameters
        ----------
        mz : ``float``
            target m/z
        mz_tol : ``float``
            m/z tolerance
        rt_bounds : ``tuple(float)``, optional
            min, max RT range to extract chromatogram

        Returns
        -------
        chrom_rts : ``np.ndarray(float)``
        chrom_ins : ``np.ndarray(float)``
            chromatogram retention time and intensity components 
        """
        rts, ins = [], []
        mz_min, mz_max = mz - mz_tol, mz + mz_tol
        for scan, srt in zip(self.ms1_scans, self.metadata.loc[self.ms1_scans, 'RetentionTime']):
            if rt_bounds is None or (srt >= rt_bounds[0] and srt <= rt_bounds[1]):  # optionally filter to only include specified RT range
                smzb, sin = np.array([self.arrays_mz.loc[scan, 'Data'], self.arrays_i.loc[scan, 'Data']])
                smz = self.mz_full[smzb.astype(np.int64)]
                rts.append(srt)
                ins.append(np.sum(sin[(smz >= mz_min) & (smz <= mz_max)]))
        return np.array([rts, ins])

    def get_msms_spectrum(self, mz, mz_tol, rt_min, rt_max, mz_bin_min, mz_bin_max, mz_bin_size):
        """
        Selects all MS2 scans with precursor m/z within tolerance of a target value and retention time 
        between specified bounds, sums spectra together returning the accumulated spectrum and some 
        metadata about how many scans were included and what precursor m/zs were included
        
        Parameters
        ----------
        mz : ``float``
            target m/z for precursor
        mz_tol : ``float``
            m/z tolerance for precursor
        rt_min : ``float``
            minimum RT for precursor
        rt_max : ``float``
            maximum RT for precursor
        mz_bin_min : ``float``
            minimum m/z for m/z binning
        mz_bin_max : ``float``
            maximum m/z for m/z binning
        mz_bin_size : ``float``
            size of bins for m/z binning
            
        Returns
        -------
        ms2_mz, ms2_i : ``np.ndarray(float)``
            mass spectrum
        n_ms2_scans : ``int``
            number of ms2 scans in spectrum
        scan_pre_mzs : ``list(float)``
            list of precursor m/zs for selected MS/MS spectrum
        """
        mz_min, mz_max = mz - mz_tol, mz + mz_tol
        # get peak precursor scans
        peak_pre_scans = self.metadata[(self.metadata['MSLevel'] == 1) & (self.metadata['RetentionTime'] >= rt_min) & (self.metadata['RetentionTime'] <= rt_max)].index.to_numpy()
        # get peak fragmentation scans
        peak_frag_scans = self.metadata[np.isin(self.metadata['PrecursorScan'], peak_pre_scans) & (self.metadata['PrecursorMonoisotopicMz'] >= mz_min) & (self.metadata['PrecursorMonoisotopicMz'] <= mz_max)].index.to_numpy()
        # precursor m/zs and count of MS2 scans
        scan_pre_mzs = self.metadata.loc[peak_frag_scans, 'PrecursorMonoisotopicMz'].tolist()
        # m/z binning parameters
        mz_bin_range = mz_bin_max - mz_bin_min
        n_bins = int((mz_bin_max - mz_bin_min) / mz_bin_size) + 1
        mz_bins = np.linspace(mz_bin_min, mz_bin_max, n_bins)
        i_bins = np.zeros(mz_bins.shape)
        mz_to_bin_idx = lambda mz: int(round((mz - mz_bin_min) / mz_bin_range * n_bins, 0))
        # accumulate MS2 scans together
        for scan in peak_frag_scans:
            for m, i in zip(self.arrays_mz.loc[scan, 'Data'], self.arrays_i.loc[scan, 'Data']):
                idx = mz_to_bin_idx(self.mz_full[m])
                if idx >= 0 and idx < n_bins:
                    i_bins[idx] += i
        # return accumulated spectrum
        return mz_bins, i_bins, len(scan_pre_mzs), scan_pre_mzs

    def get_tic(self):
        """
        TIC

        Returns
        -------
        tic_rts : ``np.ndarray(float)``
        tic_ins : ``np.ndarray(float)``
            TIC retention time and intensity components 
        """
        tic_x = self.metadata[self.metadata['PrecursorScan'] == 0].loc[:, 'RetentionTime'].to_numpy()
        tic_y = self.metadata[self.metadata['PrecursorScan'] == 0].loc[:, 'TIC'].to_numpy()
        return tic_x, tic_y

    def get_pre_mzs(self):
        """
        returns the set of unique m/z values for all MS/MS scan precursors in this file

        Returns
        -------
        pre_mzs : ``set(float)``
            sorted unique precursor m/zs
        """
        return sorted(set(self.metadata.loc[self.ms2_scans, 'PrecursorMonoisotopicMz'].tolist()))


class MsmsReaderDdaCachedMs1(MsmsReaderDda):
    """
    _MSMSReaderDDA with all arrays_mz and arrays_i data pre-read into memory to reduce 
    disk access, applicable primarily to extracting chromatograms (MS1) since that takes 
    much longer than MS2 spectra

    This takes up too much memory when there are multiple processes...

    Attributes
    ----------
    h5 : ``h5py.File``
        interface to MZA file (HDF5)
    f : ``str``
        path to MZA file
    metadata : ``pandas.DataFrame``
        metadata dataframe
    arrays_mz : ``pandas.DataFrame``
        m/z dataframe
    arrays_i : ``pandas.DataFrame``
        intensity dataframe
    ms1_scans : ``numpy.ndarray(int)``
        array of MS1 scans
    min_rt, max_rt : ``float``
        minimum/maximum retention time
    arrays_mz_cached : ``dict(int:array(float))``
    arrays_i_cached : ``dict(int:array(float))``
        m/z and intensity components of cached mass spectra (MS1 level), mapped to scan number 
    """

    def __init__(self, mza_file, drop_scans=None):
        """
        initialize the reader
    
        Parameters
        ----------
        mza_file : ``str``
            path to MZA file
        drop_scans : ``list(int)``, optional
            list of scans to drop from the file, can be None if there are not any to drop
        """
        super().__init__(mza_file, drop_scans=drop_scans)
        self.arrays_mz_cached = {}
        self.arrays_i_cached = {}
        for scan in self.ms1_scans:
            scan_mz = self.arrays_mz.loc[scan, 'Data']
            scan_i = self.arrays_i.loc[scan, 'Data']
            a_mz, a_i = np.zeros(shape=(2, scan_mz.size))
            scan_mz.read_direct(a_mz)
            scan_i.read_direct(a_i)
            self.arrays_mz_cached[scan] = self.mz_full[a_mz.astype(np.int64)]
            self.arrays_i_cached[scan] = a_i

    def close(self):
        """
        free up memory from cached data then close connection to the MZA file
        """
        del self.arrays_mz_cached
        del self.arrays_i_cached
        super().close()

    def get_chrom(self, mz, mz_tol):
        """
        Select a chromatogram (MS1 only) for a target m/z with specified tolerance

        Parameters
        ----------
        mz : ``float``
            target m/z
        mz_tol : ``float``
            m/z tolerance

        Returns
        -------
        chrom_rts : ``np.ndarray(float)``
        chrom_ins : ``np.ndarray(float)``
            chromatogram retention time and intensity components 
        """
        rts, ins = [], []
        mz_min, mz_max = mz - mz_tol, mz + mz_tol
        for scan, rt in zip(self.ms1_scans, self.metadata.loc[self.ms1_scans, 'RetentionTime']):
            smz, sin = self.arrays_mz_cached[scan], self.arrays_i_cached[scan]
            rts.append(rt)
            ins.append(np.sum(sin[(smz >= mz_min) & (smz <= mz_max)]))
        return np.array([rts, ins])
