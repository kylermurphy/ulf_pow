
# -*- coding: utf-8 -*-
"""

A set of functions to calculate the  power spectral density (PSD) from daily 
time series data from various ground-based magnetometer arrays. 


Examples
-------


Attributes
----------
local_dir : str
    Directory where the PSD Pandas dataframe is written. 
    Directory structure is local_dir/year/month/day/. 
    This is to limit the number of files in each directory for
    faster searching. 

Notes
-----

"""
import pandas as pd
import numpy as np
import scipy.signal as signal
import os

import gmag.arrays.image as image
import gmag.arrays.carisma as carisma
import gmag.utils as g_utils

#TODO: add config parser for data dir
#TODO: add additional function to the different power functions? 
# By station, entire array? array_pow(...)?

local_dir = "E:\\data\\magnetometer\\psd\\"

def image_pow(site: str = ['AND'], sdate='2010-01-01', edate=None, ndays=1):
    """Calculate power for magnetometer stations
    from the IMAGE array. 

    Hourly Power Spectral Density (PSD) is calculated
    for the D component of the ground-based magnetoemeter.

    Power is stored in a Panda's dataframe and written
    out as compressed csv file.
    
    Parameters
    ----------
    site : str, optional
        Magnetometer stations to be processed, by default ['AND']
    sdate : str, optional
        Start date to process power, by default '2010-01-01'
    edate : [type], optional
        End date for calculating power, by default None
    ndays : int, optional
        Number of days to process following start date, by default 1.
        This is ignored if edate is specified.
    """
    if type(site) is str:
        site = [site]

    # get a list of days to loop over
    if edate is not None:
        d_arr = pd.Series(pd.date_range(start=sdate, end=edate, freq='D'))
    else:
        d_arr = pd.Series(pd.date_range(start=sdate, periods=ndays, freq='D'))

    for di, dt in d_arr.iteritems():
        # load image data
        i_dat = image.load(site=site, sdate=dt, ndays=1)
        if i_dat is None:
            continue
        # get nominal resolution
        # the diff that occurs the most
        res = (pd.Series(i_dat.index[1:]) -
               pd.Series(i_dat.index[:-1])).value_counts()
        res = res.index[0].total_seconds()

        for ss in site:
            try:
                d = i_dat[ss.upper()+'_D'].copy()
            except KeyError:
                print('East/West data not returned for {0}'.format(ss))
                continue

            print('Calculating power for {0}'.format(ss))
            hour_pow(d=d,
                     sdate=dt,
                     decl=i_dat.iloc[0][ss.upper()+'_declination'],
                     cgm_lat=i_dat.iloc[0][ss.upper()+'_cgmlat'],
                     cgm_lon=i_dat.iloc[0][ss.upper()+'_cgmlon'],
                     l_shell=i_dat.iloc[0][ss.upper()+'_lshell'],
                     site=ss.upper(),
                     res=res)


def carisma_pow(site: str = ['GILL'], sdate='2010-01-01', edate=None, ndays=1):
    """Calculate power for magnetometer stations
    from the CARISMA array. 

    Hourly Power Spectral Density (PSD) is calculated
    for the D component of the ground-based magnetoemeter. 

    Power is stored in a Panda's dataframe and written
    out as compressed csv file.
    
    Parameters
    ----------
    site : str, optional
        Magnetometer stations to be processed, by default ['GILL']
    sdate : str, optional
        Start date to process power, by default '2010-01-01'
    edate : [type], optional
        End date for calculating power, by default None
    ndays : int, optional
        Number of days to process following start date, by default 1.
        This is ignored if edate is specified.
    """
    if type(site) is str:
        site = [site]

    # get a list of days to loop over
    if edate is not None:
        d_arr = pd.Series(pd.date_range(start=sdate, end=edate, freq='D'))
    else:
        d_arr = pd.Series(pd.date_range(start=sdate, periods=ndays, freq='D'))

    for di, dt in d_arr.iteritems():
        for ss in site:
            # load carisma data
            i_dat = carisma.load(site=ss.upper(), sdate=dt, ndays=1)

            if i_dat is None:
                continue
            try:
                d = i_dat[ss.upper()+'_D'].copy()
            except KeyError:
                print('East/West data not returned for {0}'.format(ss))
                continue
            # get nominal resolution
            res = (pd.Series(i_dat.index[1:]) -
                   pd.Series(i_dat.index[:-1])).value_counts()
            res = res.index[0].total_seconds()
            # calculate power
            hour_pow(d=d,
                     sdate=dt,
                     decl=i_dat.iloc[0][ss.upper()+'_declination'],
                     cgm_lat=i_dat.iloc[0][ss.upper()+'_cgmlat'],
                     cgm_lon=i_dat.iloc[0][ss.upper()+'_cgmlon'],
                     l_shell=i_dat.iloc[0][ss.upper()+'_lshell'],
                     site=ss.upper(),
                     res=res)


def hour_pow(d, sdate, decl, cgm_lat, cgm_lon, l_shell, site, res, flen=3600, s_thresh=25.):
    """Calculate hourly PSD from a time series passed by
    one of the array_pow functions. 

    PSD is calculated using the Numpy singal processing
    function signal.periodogram.

    PSD is calculated assuming a full day of data is passed to 
    the function. 24 hourly PSDs are produced.

    PSD is saved in a Pandas dataframe and and written out
    into a compressed csv. 

    Parameters
    ----------
    d : [type]
        Daily magnetometer time series to calculate hourly power
        for
    sdate : [type]
        Date of time series for creating directory structure and
        filenames
    decl : [type]
        Station declination.
    cgm_lat : [type]
        Station corrected geomagnetic latitude.
    cgm_lon : [type]
        Station corrected geomagnetic longitude.
    l_shell : [type]
        Station L-shell.
    site : [type]
        Station code.
    res : [type]
        Time series resolution, required for calculating frequency
        resolution, and calculating the window length for each 
        time series
    flen : int, optional
        Assumed 1 hour time series, by default 3600. 
        Variation of this will allow for PSD to be calculated for
        different windows.
    s_thresh : [type], optional
        A threshold for identify spikes in the time series data, by default 25.
        Simply compares tn and tn-1. If the difference is greater then the threshold
        replace with NaN. 
        This will produce a NaN PSD array when calculating power.
    """
    flen = int(3600 / res)
    fs = 1. / res

    # dummy spectra for size
    fp, p_spec = signal.periodogram(np.arange(
        0, flen), fs, scaling='density', detrend='constant', window='hann')

    # dummy/nan power
    p_nan = p_spec[:]
    p_nan[:] = np.nan

    # check for spikes
    d_diff = d - d.shift()
    d[abs(d_diff) > s_thresh] = np.nan
    if any(abs(d_diff) > s_thresh):
        print("Spikes found: {0}".format(sdate))
    # create a data frame to store
    #  the power spectral density data
    pf = pd.DataFrame({'freq': fp,
                       'station': site,
                       'decl': float(decl),
                       'L': float(l_shell),
                       'cgm_lat': float(cgm_lat),
                       'cgm_lon': float(cgm_lon)})

    # calculate power spectral density
    for j in np.arange(0, 24):
        d_fft = d[j * flen:(j + 1) * flen]
        if any(d_fft.isnull()):
            print("NaN found at hour {0}".format(j))
            p_spec = p_nan[:]
        else:
            fp, p_spec = signal.periodogram(
                d_fft, fs, scaling='density', detrend='constant', window='hann')
        # fill hour power
        pf['H{0:02d}'.format(j)] = p_spec

    # create save directory and
    # check directory exists
    f_d = os.path.join(local_dir, '{0:04d}'.format(
        sdate.year), '{0:02d}'.format(sdate.month), '{0:02d}'.format(sdate.day))
    if not os.path.exists(f_d):
        os.makedirs(f_d)

    # print power to csv
    fn = os.path.join(f_d, '{0:04d}{1:02d}{2:02d}{3}_psd.txt.gz'.format(
        sdate.year, sdate.month, sdate.day, site))
    pf.to_csv(fn, index=False, float_format="%E",
              na_rep='NaN', compression='gzip')


