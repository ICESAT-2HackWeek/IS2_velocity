#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np

def fill_seg_ids(data,dx=20):
    r"""Fill the along-track vector so that there are no skipped points

    Parameters
    ------
    data : dict
        data dictionary
    dx : float; default 20
         step between points (20 m is the default for atl06)

    Output
    ------
    data : dict
    """

    variables = ['x_atc_full','h_full','lats_full','lons_full','x_full','y_full']
    for var in variables:
        data[var] = {}

    # make a monotonically increasing x vector
    for cycle in data['x_atc'].keys():
        for var in variables:
            data[var][cycle] = {}
        for beam in data['x_atc'][cycle].keys():
            # indices starting at zero, using the segment_id field, so any skipped segment will be kept in correct location
            seg_ids = data['segment_ids'][cycle][beam]
            ind = seg_ids - np.nanmin(seg_ids)

            data['x_atc_full'][cycle][beam] = np.arange(np.max(ind)+1) * dx + data['x_atc'][cycle][beam][0]
            data['h_full'][cycle][beam] = np.zeros(np.max(ind)+1) + np.NaN
            data['h_full'][cycle][beam][ind] = data['h_li'][cycle][beam]

            data['lats_full'][cycle][beam] = np.zeros(np.shape(data['x_atc_full'][cycle][beam])) * np.nan
            data['lats_full'][cycle][beam][ind] = data['latitude'][cycle][beam]
            data['lons_full'][cycle][beam] = np.zeros(np.shape(data['x_atc_full'][cycle][beam])) * np.nan
            data['lons_full'][cycle][beam][ind] = data['longitude'][cycle][beam]
            data['x_full'][cycle][beam] = np.zeros(np.shape(data['x_atc_full'][cycle][beam])) * np.nan
            data['x_full'][cycle][beam][ind] = data['x'][cycle][beam]
            data['y_full'][cycle][beam] = np.zeros(np.shape(data['x_atc_full'][cycle][beam])) * np.nan
            data['y_full'][cycle][beam][ind] = data['y'][cycle][beam]

            ## save the segment id's themselves, with gaps filled in
            data['segment_ids'][cycle][beam] = np.zeros(np.max(ind)+1) + np.NaN
            data['segment_ids'][cycle][beam][ind] = seg_ids

    return data

# -------------------------------------------------------------------------------------------

def interpolate_nans(data):
    r"""Linear interpolation between nan values

    Parameters
    ------
    data : dict
        data dictionary

    Output
    ------
    data : dict
    """

    Warning("Be careful interpolating, linear interpolation could lead to misleading correlations.")

    data['h_full_interp'] = {}
    for cycle in data['h_full'].keys():
        data['h_full_interp'][cycle] = {}
        for beam in data['h_full'][cycle].keys():
            ### interpolate nans in pandas
            interp_hold = {'x_atc_full': data['x_atc_full'][cycle][beam], 'h_full': data['h_full'][cycle][beam]}
            df = pd.DataFrame(interp_hold, columns = ['x_atc_full','h_full'])
            # linear interpolation for now
            df['h_full'].interpolate(method = 'linear', inplace = True)
            data['h_full_interp'][cycle][beam] = df['h_full'].values
    return data

# -------------------------------------------------------------------------------------------

def filt(data, dx=None, filter_type = 'running_average', running_avg_window = None):
    r"""Smoothing filter to remove high-frequency noise

    Parameters
    ------
    data : dict
        data dictionary
    dx : float; default None
        spatial step between points; if None use median
    filter_type : str; default running_average
        choice on how to filter the data
    running_avg_window : int
        window size of the filter (number of samples)

    Output
    ------
    data : dict
    """

    data['h_filt'] = {}
    for cycle in data['h_full_interp'].keys():
        data['h_filt'][cycle] = {}
        for beam in data['h_full_interp'][cycle].keys():
            if filter_type == 'running_average': #
                if not running_avg_window: # if no window size was given, run a default 100 m running average smoother
                    running_avg_window = 100
                if dx is None:
                    dx = np.nanmedian(data['x_atc_full'][cycle][beam][1:] - data['x_atc_full'][cycle][beam][:-1])
                filt = np.ones(int(np.round(running_avg_window / dx)))
                data['h_filt'][cycle][beam] = (1/len(filt)) * np.convolve(filt, data['h_full_interp'][cycle][beam], mode = 'same')

            #TODO: add bandpass, highpass
            else:
                data['h_filt'][cycle][beam] = data['h_full_interp'][cycle][beam]

    return data

# -------------------------------------------------------------------------------------------

def differentiate(data,dx=None):
    r"""Slope calculation

    Parameters
    ------
    data : dict
        data dictionary
    dx : float; default None
        spatial step between points; if None use median

    Output
    ------
    data : dict
    """

    data['h_li_diff'] = {}
    for cycle in data['h_filt'].keys():
        data['h_li_diff'][cycle] = {}
        for beam in data['h_filt'][cycle].keys():
            if dx is None:
                dx = np.nanmedian(data['x_atc_full'][cycle][beam][1:] - data['x_atc_full'][cycle][beam][:-1])
            data['h_li_diff'][cycle][beam] = np.gradient(data['h_filt'][cycle][beam], dx)
    return data

