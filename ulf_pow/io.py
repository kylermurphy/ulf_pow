# -*- coding: utf-8 -*-
"""

This module supports some basic reading and writing of datafiles produced by pow.py.
It's meant to simplify reading in of the text files produced in pow.py and create
a simpler file format using hdf5 files and PyTables which allow for rapid reading as
well as searching data columns when reading in the data. For example reading in specific 
chuncks, date ranges, l-shells, stations, omni set; this can signficantly save memory
analysing the data.

get_spec reads in the text files for a fixed frequency range.

sum_psd reads in and sums/integrates the text files for a fixed frequency range.

ulf_database calls either sum_psd or get_spec for a date range and number of stations
and outputs the data into an hdf5 file. This is meant to save time so that text files
don't need to be continously read in for analysis. 


Attributes
----------

local_dir : str
    Directory where the PSD files are stored.

"""

import pandas as pd
import numpy as np
import scipy.signal as signal
import os

from datetime import timedelta

from gmag.omni import omni
from gmag import utils

local_dir = "D:\\data\\magnetometer\\psd\\"

def get_spec(site, sdate,  edate=None, ndays=1,
                psd_dir=local_dir,
                fmin=0.80, fmax=15.01, add_omni: bool=False, omni_lags=[0],
                verbose: bool=False):
    """Returns a Pandas DataFrame containing the PSD spectrum, station name, 
    and magnetic local time (MLT) of the station for a set of dates. 

    Will loop through multiple stations if provided.
    
    Parameters
    ----------
    site : str
        String or array of strings with the station names 
        which the spectra will be loaded
    sdate : str of datetime-like
        Initial day to be loaded
    edate : [str of datetime-like
        Final day to be loaded
    ndays : int, optional
        Number of days to load if edate is not defined, by default 1
    psd_dir : str, optional
        Directory of PSD files, by default local_dir
    fmin : float, optional
        Minimum frequency (mHz) to load, by default 0.80
    fmax : float, optional
        Maximum frequency (mHz) to load, by default 15.01
    add_omni : bool, optional
        Add hourly omni data to data frame, by default False
    omni_lags : list, optional
        Add lagged omni data; list is a list of hourly lags to include,
        by default [0], only include the current hours omni data.
        This is attached to each spectra or row in the data frame
        and so the returned DataFrame can get large quite quickly if 
        looking at large lag times.
    verbose : bool, optional
        Print some simple messages, by default False
    
    Returns
    -------
    Pandas DataFrame
        DataFrame containing the spectrum, station, MLT, frequency resolution
    """
    if add_omni:
        if type(omni_lags) is int:
            omni_lags = [omni_lags]
        elif  type(omni_lags) is np.ndarray:
            omni_lags = omni_lags.tolist()
        
        max_l = int(max(omni_lags)+1)
        o_dat = omni.load(pd.to_datetime(sdate)-timedelta(hours=max_l), 
                          pd.to_datetime(edate)+timedelta(hours=max_l))
        
        print('Adding omni data with lags:{0}'.format(omni_lags))

    if type(site) is str:
        site = [site.upper()]

    # get a list of days to loop over
    if edate is not None:
        d_arr = pd.Series(pd.date_range(start=sdate, end=edate, freq='D'))
    else:
        d_arr = pd.Series(pd.date_range(start=sdate, periods=ndays, freq='D'))


    # convert fmin and fmax to Hz
    fmin = fmin/1000.
    fmax = fmax/1000.


    df_r = pd.DataFrame() 

    for stn in site:
        yr = -1
        for di, dt in d_arr.iteritems(): 
            df_t = pd.Series(pd.date_range(dt, freq='H', periods=24))

            fn = '{0:04d}{1:02d}{2:02d}{3}_psd.txt.gz'.format(
                dt.year, dt.month, dt.day, stn.upper())

            fn = os.path.join(psd_dir,
                              '{0:04d}'.format(dt.year),
                              '{0:02d}'.format(dt.month),
                              '{0:02d}'.format(dt.day),
                              fn)

            # read in file, create an empty
            # data frame if file is not found
            try:
                df_in = pd.read_csv(fn, compression='gzip', header=0)
            except FileNotFoundError:
                if verbose:
                    print("File not found {0}".format(fn))
                continue

            if verbose:
                print("File loaded {0}".format(fn))
            # drop unnecessary columns
            df_in = df_in[(df_in.freq >= fmin) & (df_in.freq <= fmax)]
            freq = df_in['freq'].copy()*1000.
            df_in = df_in.drop(labels=['freq','station','decl','L','cgm_lat','cgm_lon'],axis=1)
            # get frequency resolution
            f_res = freq.iloc[1]-freq.iloc[0]
            # trasnpose array and drop new index column
            # convert from nT^2/Hz  to nT^2/mHz
            df_in = df_in.transpose().reset_index()
            df_in = df_in.drop(labels=['index'],axis=1)
            df_in = df_in.div(1000.)
            # add time to the array
            df_in['t'] = df_t
            df_in['stn'] = stn.upper()
            df_in['f res mHz'] = f_res
            df_in = df_in.set_index('t')
            # add MLT to the array
            if dt.year != yr:
                print(stn, dt.year)
                stn_c = utils.load_station_coor(param=stn, year=dt.year)
            mlt = np.arange(0, 24)
            mlt = (mlt+24-float(stn_c['mlt_midnight'])) % 24.
            df_in['mlt'] = mlt

            df_r = df_r.append(df_in)
            yr = dt.year

        # append read in data
        if add_omni:
            df_in = df_in.join(o_dat, how='left')

        

    return df_r




def sum_psd(site, sdate, edate=None, ndays=1,
                psd_dir=local_dir,
                fmin=0.80, fmax=15.01, add_omni: bool=False, omni_lags=[0],
                verbose: bool=False):
    """Returns a Pandas DataFrame containing the integrated and summed PSD.

    Will loop through multiple stations if provided.
    
    
    Parameters
    ----------
    site : str
        String or array of strings with the station names 
        which the spectra will be loaded
    sdate : str of datetime-like
        Initial day to be loaded
    edate : [str of datetime-like
        Final day to be loaded
    ndays : int, optional
        Number of days to load if edate is not defined, by default 1
    psd_dir : str, optional
        Directory of PSD files, by default local_dir
    fmin : float, optional
        Minimum frequency (mHz) to load, by default 0.80
    fmax : float, optional
        Maximum frequency (mHz) to load, by default 15.01
    add_omni : bool, optional
        Add hourly omni data to data frame, by default False
    omni_lags : list, optional
        Add lagged omni data; list is a list of hourly lags to include,
        by default [0], only include the current hours omni data.
        This is attached to each spectra or row in the data frame
        and so the returned DataFrame can get large quite quickly if 
        looking at large lag times.
    verbose : bool, optional
        Print some simple messages, by default False
    
    Returns
    -------
    Pandas DataFrame
        Date frame with integrated and summed PSD between fmin and fmax,
        MLT, L-shell, station.
    """
    if add_omni:
        if type(omni_lags) is int:
            omni_lags = [omni_lags]
        elif  type(omni_lags) is np.ndarray:
            omni_lags = omni_lags.tolist()
        
        max_l = int(max(omni_lags)+1)
        o_dat = omni.load(pd.to_datetime(sdate)-timedelta(hours=max_l), 
                          pd.to_datetime(edate)+timedelta(hours=max_l))
        
        print('Adding omni data with lags:{0}'.format(omni_lags))

    if type(site) is str:
        site = [site.upper()]

        # get a list of days to loop over
    if edate is not None:
        d_arr = pd.Series(pd.date_range(start=sdate, end=edate, freq='D'))
    else:
        d_arr = pd.Series(pd.date_range(start=sdate, periods=ndays, freq='D'))

    df_r = pd.DataFrame()

    for stn in site:
        df_psd = pd.DataFrame()
        yr = -1
        for di, dt in d_arr.iteritems():
            # hourly array for time axis
            df_t = pd.Series(pd.date_range(dt, freq='H', periods=24))

            fn = '{0:04d}{1:02d}{2:02d}{3}_psd.txt.gz'.format(
                dt.year, dt.month, dt.day, stn.upper())

            fn = os.path.join(psd_dir,
                              '{0:04d}'.format(dt.year),
                              '{0:02d}'.format(dt.month),
                              '{0:02d}'.format(dt.day),
                              fn)

            # read in file, create an empty
            # data frame if file is not found
            try:
                df_in = pd.read_csv(fn, compression='gzip', header=0)
            except FileNotFoundError:
                if verbose:
                    print("File not found: {0}".format(fn))
                sp = pd.DataFrame({'t': df_t})
                sp['ipsd'] = np.nan
                sp['spsd'] = np.nan
                sp['mlt'] = np.nan
                sp['lshell'] = np.nan

                df_psd = df_psd.append(sp, ignore_index=True, sort=True)
                continue

            if verbose:
                print("File loaded: {0}".format(fn))
            # read in station coordinates
            if dt.year != yr:
                print(stn, dt.year)
                stn_c = utils.load_station_coor(param=stn, year=dt.year)

            # convert frequency to mHz
            df_in['freq'] = df_in['freq'] * 1000.

            # find frequency range
            gf = df_in['freq'].between(fmin, fmax)
            gd = df_in[gf]

            # sum power and convert from
            # nT^2/Hz to nT^2/mHz
            sx = gd.sum(axis=0, min_count=gd.shape[0])['H00':'H23']
            sx = sx / 1000.
            # integrate power
            ix = sx / (gd['freq'].max() - gd['freq'].min())

            # fill data frame
            mlt = np.arange(0, 24)
            mlt = (mlt+24-float(stn_c['mlt_midnight'])) % 24.
            sp = pd.DataFrame({'t': df_t,
                               'ipsd': ix.values,
                               'spsd': sx.values,
                               'lshell': float(stn_c['lshell']),
                               'mlt':  mlt})

            df_psd = df_psd.append(sp, ignore_index=True, sort=True)
            yr = dt.year

        if df_psd.empty:
            continue

        df_psd['ipsd'] = pd.to_numeric(df_psd['ipsd'])
        df_psd['spsd'] = pd.to_numeric(df_psd['spsd'])
        df_psd['stn'] = stn.upper()
        df_psd = df_psd.sort_values(by=['t']).reset_index(drop=True)
        df_psd = df_psd.set_index('t')
        df_psd = df_psd.dropna(subset=['ipsd'])

        # add omni data and lagged
        # omni if needed to each station
        if add_omni:
            o_t = o_dat.index
            o_c = o_dat.columns
            for lags in omni_lags:
                if lags != 0:
                    o_dat['t'] = o_t + timedelta(hours=lags)
                    o_dat = o_dat.set_index('t')
                    o_dat.columns = o_c+'_{0}'.format(lags)
                df_psd = df_psd.join(o_dat, how='left', sort=True)

        # getting pandas warning here
        # do we need to sort here? 
        if df_r.empty:
            df_r = df_psd.copy()
        else:
            df_r = df_r.append(df_psd)

    return df_r



def ulf_database(site=None, sdate='2012-01-01', edate=None, ndays=1,
                psd_dir=local_dir,
                db_type = 'sum',
                fmin=0.80, fmax=15.01, add_omni: bool=False, omni_lags=[0],
                verbose: bool=False, 
                psd_f: str=None):
    """Simple function to load and output data from sum_psd and get_spec into
    hdf5 files. 

    Purpose was to simplify sorting through large amounts of data (dates range, 
    number of stations) and make loading data for following analysis simpler. 

    Utility can be easily replicated by just calling sum_psd or get_spec. 
    However, the function here does not pass an array of stations to either
    function, rather it loops through each station and appends the data to 
    the saved hdf5 file. This can save memory when creating a database of a 
    large number of stations or an extended date range.
    
    Parameters
    ----------    

    site : str
        String or array of strings with the station names 
        which the spectra will be loaded
    sdate : str of datetime-like
        Initial day to be loaded
    edate : str of datetime-like
        Final day to be loaded
    ndays : int, optional
        Number of days to load if edate is not defined, by default 1
    psd_dir : str, optional
        Directory of PSD files, by default local_dir
    db_type : str, optional
        Type of data to get, by default 'sum', calls sum_psd otherwise
        will call get_spec
    fmin : float, optional
        Minimum frequency (mHz) to load, by default 0.80
    fmax : float, optional
        Maximum frequency (mHz) to load, by default 15.01
    add_omni : bool, optional
        Add hourly omni data to data frame, by default False
    omni_lags : list, optional
        Add lagged omni data; list is a list of hourly lags to include,
        by default [0], only include the current hours omni data.
        This is attached to each spectra or row in the data frame
        and so the returned DataFrame can get large quite quickly if 
        looking at large lag times.
    verbose : bool, optional
        Print some simple messages, by default False
    psd_f : str, optional
        File name to save.
    
    Returns
    -------
    [type]
        [description]
    """
    # a set of ground based site currently 
    # begin used in larger analysis
    # sites are from CARISMA and IMAGE
    if site is None:
        site = ['BACK', 'DAWS', 'FCHP', 'FSIM', 'FSMI',
                'GILL', 'GULL', 'ISLL', 'LGRR', 'MCMU',
                'MSTK', 'OSAK', 'OXFO', 'PINA', 'POLS',
                'RABB', 'THRF', 'VULC', 'WEYB', 'WGRY',
                'ABK', 'ALT', 'AND', 'DOB', 'DON',
                'HAN', 'HAR', 'IVA', 'JCK', 'KAR',
                'KAU', 'KEV', 'KIL', 'KIR', 'LEK',
                'LOZ', 'LYC', 'MAS', 'MEK', 'MUO',
                'NUR', 'OUJ', 'PEL', 'RAN', 'RST',
                'RVK', 'SOD', 'SOL', 'SOR', 'TAR',
                'TRO', 'UPS'] 

    # determine what type of file is going
    # to be out put
    # either summed psd
    # or the power spectra
    # and set filename if not provided
    if db_type.lower() == 'sum':
        ulf_fun = sum_psd
        dat_col = True
        db_type = 'sum_psd'
        if psd_f is None:
            psd_f = 'ulf_psd_summed.h5'
    else:
        ulf_fun = get_spec
        dat_col = ['t','stn','mlt']
        db_type = 'get_spec'
        if psd_f is None:
            psd_f = 'ulf_psd_spec.h5'
    
    # loop over stations and grab PSD
    pow_dat = pd.DataFrame()

    for stn in site:
        print("Returning {0}:".format(db_type))
        print(site)

        pow_dat = ulf_fun(stn, sdate, edate=edate, ndays=ndays, 
            psd_dir=psd_dir, 
            fmin=0.80, fmax=15.01, add_omni=add_omni, omni_lags=omni_lags,
            verbose=verbose)

        # hdf5 does not like NaN, this causes problem when querying the data
        # this should fix that
        pow_dat = pow_dat.fillna(-999999.999)
        # write out data
        pow_dat.to_hdf(os.path.join(psd_dir, psd_f), key='pow_dat',
                    mode='a', append=True, 
                    format='table', data_columns=dat_col,                  
                    complevel=9, complib='blosc:snappy')

    return pow_dat








