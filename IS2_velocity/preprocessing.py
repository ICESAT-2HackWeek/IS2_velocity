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

    # make a monotonically increasing x vector
    for cycle in data['x_atc'].keys():
        for beam in data['x_atc'][cycle].keys():
            # indices starting at zero, using the segment_id field, so any skipped segment will be kept in correct location
            seg_ids = data['segment_ids'][cycle][beam]
            ind = seg_ids - np.nanmin(seg_ids)

            data['x_atc'][cycle][beam] = np.arange(np.max(ind)+1) * dx + data['x_atc'][cycle][beam][0]

            h_li_hold = np.zeros(np.max(ind)+1) + np.NaN
            h_li_hold[ind] = data['h_li'][cycle][beam]
            data['h_li'][cycle][beam] = h_li_hold

            lats_hold = np.zeros(np.shape(data['x_atc'][cycle][beam])) * np.nan
            lats_hold[ind] = data['latitude'][cycle][beam]
            data['latitude'][cycle][beam] = lats_hold

            lons_hold = np.zeros(np.shape(data['x_atc'][cycle][beam])) * np.nan
            lons_hold[ind] = data['longitude'][cycle][beam]
            data['longitude'][cycle][beam] = lons_hold

            xs_hold = np.zeros(np.shape(data['x_atc'][cycle][beam])) * np.nan
            xs_hold[ind] = data['x'][cycle][beam]
            data['x'][cycle][beam] = xs_hold

            ys_hold = np.zeros(np.shape(data['x_atc'][cycle][beam])) * np.nan
            ys_hold[ind] = data['y'][cycle][beam]
            data['y'][cycle][beam] = ys_hold

            ## save the segment id's themselves, with gaps filled in
            data['segment_ids'][cycle][beam] = np.zeros(np.max(ind)+1) + np.NaN
            data['segment_ids'][cycle][beam][ind] = seg_ids

    data['flags'].append('fill_seg_ids')

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

    # to save raw data
    data['h_li_RAW'] = {}

    for cycle in data['h_li'].keys():
        # to save raw data
        data['h_li_RAW'][cycle] = {}

        for beam in data['h_li'][cycle].keys():
            # to save raw data
            data['h_li_RAW'][cycle][beam] = data['h_li'][cycle][beam]

            ### interpolate nans in pandas
            interp_hold = {'x_atc': data['x_atc'][cycle][beam], 'h_li': data['h_li'][cycle][beam]}
            df = pd.DataFrame(interp_hold, columns = ['x_atc','h_li'])
            # linear interpolation for now; forward and backward, so nans at start and end are filled also
            df['h_li'].interpolate(method = 'linear', limit_direction = 'forward', inplace = True)
            df['h_li'].interpolate(method = 'linear', limit_direction = 'backward', inplace = True)
            data['h_li'][cycle][beam] = df['h_li'].values

    data['flags'].append('interp_nans')
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

    for cycle in data['h_li'].keys():
        for beam in data['h_li'][cycle].keys():
            if filter_type == 'running_average': #
                if not running_avg_window: # if no window size was given, run a default 100 m running average smoother
                    running_avg_window = 100
                if dx is None:
                    dx = np.nanmedian(data['x_atc'][cycle][beam][1:] - data['x_atc'][cycle][beam][:-1])
                filt = np.ones(int(np.round(running_avg_window / dx)))
                data['h_li'][cycle][beam] = (1/len(filt)) * np.convolve(filt, data['h_li'][cycle][beam], mode = 'same')

            #TODO: add options bandpass, highpass
            else:
                data['h_li'][cycle][beam] = data['h_li'][cycle][beam]

    data['flags'].append('filt_'+filter_type)

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
    for cycle in data['h_li'].keys():
        data['h_li_diff'][cycle] = {}
        for beam in data['h_li'][cycle].keys():
            if dx is None:
                dx = np.nanmedian(data['x_atc'][cycle][beam][1:] - data['x_atc'][cycle][beam][:-1])
            data['h_li_diff'][cycle][beam] = np.gradient(data['h_li'][cycle][beam], dx)

    data['flags'].append('differentiate')

    return data

