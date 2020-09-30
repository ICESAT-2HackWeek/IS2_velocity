#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy.signal import correlate
from astropy.time import Time

# -------------------------------------------------------------------------------------------

def velocity(x_in,dh1,dh2,dt,output_xs,search_width=100,segment_length=2000,dx=20,corr_threshold=.65):
    r"""Calculate along-track velocity by correlating the along-track elevation between repeat acquisitions.

    Parameters
    ------
    x_in : array
           along-track distance
    dh1 : array
          surface slope at cycle 1
    dh2 : array
          surface slope at cycle 2
    dt : float
         time difference between cycles
    output_xs : array
                along-track distances at which the velocities will be calculated and output
    search_width : float
                   the maximum distance which the stencil array will be moved to look for a good correlation
    segment_length : float
                     the length of the array to be correlated
    dx : float
         spacing between points in the x array
    corr_threshold : float
                     minimum correlation to be considered an acceptable fit

    Output
    ------
    velocities : array
                 calculated velocities corresponding to locations output_xs
    correlations : array
                   calculated correlations for each of the velocity points
    """

    # create output arrays to be filled
    velocities = np.nan*np.ones_like(output_xs)
    correlations = np.nan*np.ones_like(output_xs)

    # iterate over all the output locations
    for xi,x_start in enumerate(output_xs):
        # find indices for the stencil
        idx1 = np.argmin(abs(x_in-x_start))
        idx2 = idx1 + int(np.round(segment_length/dx))
        # cut out a wider chunk of data at time t2 (second cycle)
        idx3 = idx1 - int(np.round(search_width/dx)) # offset on earlier end by # indices in search_width
        idx4 = idx2 + int(np.round(search_width/dx)) # offset on later end by # indices in search_width

        # get the two arrays that will be correlated to one another
        dh_stencil = dh1[idx1:idx2]
        dh_wide = dh2[idx3:idx4]

        # correlate old with newer data
        corr = correlate(dh_stencil, dh_wide, mode = 'valid', method = 'direct')
        norm_val = np.sqrt(np.sum(dh_stencil**2)*np.sum(dh_wide**2)) # normalize so values range between 0 and 1
        corr = corr / norm_val
        lagvec = np.arange(- int(np.round(search_width/dx)), int(search_width/dx) +1,1)# for mode = 'valid'
        shift_vec = lagvec * dx

        if all(np.isnan(corr)):
            continue
        else:
            correlation_value = np.nanmax(corr)
            if correlation_value >= corr_threshold:
                idx_peak = np.arange(len(corr))[corr == correlation_value][0]
                best_shift = shift_vec[idx_peak]
                velocities[xi] = best_shift/(dt/365)
                correlations[xi] = correlation_value
            else:
                velocities[xi] = np.nan
                correlations[xi] = correlation_value

    return velocities,correlations

# -------------------------------------------------------------------------------------------

def smooth_and_diff(x_in,h_in,win=None):
    r"""Smooth and differentiate an along-track height array.

    Parameters
    ------
    x_in : array
           along-track distance
    h_in : array
           along-track elevation
    win : float; default None
          smoothing window (in meters)
          for example, if win=60 m, the smoothing window is a 3 point running average smoothed dataset, because each point is 20 m apart

    Output
    ------
    h : array
        along-track height (smoothed if win is not None)
    dh : array
         surface slope
    """

    # Default is no smoothing
    if win is None:
        # calculate the surface slope
        dh = np.gradient(h_in,x_in)
        return h_in,dh

    # running average smoother / filter
    else:
        # calculate the step between x points
        dx = np.mean(np.gradient(x_in))
        # number of points (meters / dx)
        N = int(np.round(win/dx))
        # smooth
        h_smooth = np.nan*np.ones_like(x_in)
        h_smooth[N//2:-N//2+1] = np.convolve(h_in,np.ones(N)/N, mode='valid')
        # calculate the surface slope
        dh = np.gradient(h_smooth,x_in)
        return h_smooth,dh

# -------------------------------------------------------------------------------------------

def fill_seg_ids(x_in,h_in,seg_ids,dx=20):
    r"""Fill the along-track vector so that there are no skipped points

    Parameters
    ------
    x_in : array
           along-track distance
    h_in : array
           along-track elevation
    seg_ids : array
              segment ids as imported with atl06_to_dict function
    dx : float; default 20
         step between points (20 m is the default for atl06)

    Output
    ------
    x_full : array
             filled array of along-track distance
    h_full : array
             filled array of alont-track elevation
    """

    # make a monotonically increasing x vector
    ind = seg_ids - np.nanmin(seg_ids) # indices starting at zero, using the segment_id field, so any skipped segment will be kept in correct location
    x_full = np.arange(np.max(ind)+1) * dx + x_in[0]
    h_full = np.zeros(np.max(ind)+1) + np.NaN
    h_full[ind] = h_in

    return x_full,h_full

# -------------------------------------------------------------------------------------------

def time_diff(D1,D2):
    r"""Get the time difference between two cycles

    Parameters
    ------
    D1 : Dictionary 1
    D2 : Dictionary 2

    Output
    ------
    dt : float
         time step
    """

    # get time strings
    t1_string = D1['data_start_utc'].astype(str) #figure out later if just picking the first one it ok
    t2_string = D2['data_start_utc'].astype(str) #figure out later if just picking hte first one it ok
    # use astropy to get time from strings
    t1 = Time(t1_string)
    t2 = Time(t2_string)
    # time step to julien days
    dt = (t2 - t1).jd

    return dt

# -------------------------------------------------------------------------------------------
