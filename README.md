# ULF Power Module

This code calculates hourly power spectral density (PSD) from ground-based magnetometer data. It uses the module ```gmag``` to read in ground-based magnetometer data and ```scipy``` to calculate PSD. The power spectra density is saved as csv files using Pandas. The code also has some basic reading and writing of the saved csv files.

- ```sum_psd``` will read in the csv PSD and calculate summed and integrated over a fixed frequency range for a given date range or number of days.
- ```get_spec``` will read in the csv PSD spectra over a fixed frequency range for a given date range or number of days.
- ```ulf_database``` will read in a blocks of PSD (stations and dates) and output an hdf5 of summed (and integrated) or spectra PSD. This simplifies analysis of the dataset. 

Each of the reading in codes can also added hourly lagged omni to the PSD. This requires the ```gmag``` and ```heliopy``` packages. 

