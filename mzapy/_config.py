"""
mzapy/_config.py

Dylan Ross (dylan.ross@pnnl.gov)

    internal module with configuration parameters
"""


# versions of the underlying .mza format that are supported
_MZA_VERSIONS_SUPPORTED = [
    'old',
    'new',
]


# all metadata headers
_ALL_METADATA_HEADERS = [
    'Scan', 'MSLevel', 'Polarity', 'Fragmentation', 'RetentionTime', 'PrecursorScan', 'PrecursorMonoisotopicMz', 
    'PrecursorCharge', 'IsolationWindowTargetMz', 'IsolationWindowWidth', 'IsolationWindowOffset', 'TIC', 'ScanLabel', 
    'IonMobilityFrame', 'IonMobilityBin', 'IonMobilityTime'
]


# specify metadata headers to cache in memory
_CACHE_METADATA_HEADERS = [
    'Scan', 'MSLevel', 'RetentionTime', 'IonMobilityFrame', 'IonMobilityBin', 'IonMobilityTime'
]
