#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy.signal import correlate
from astropy.time import Time
import os,re,h5py,pyproj
from IS2_velocity.extract_alongtrack import get_measures_along_track_velocity

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

def filt(x1, h1, dx, filter_type = 'running_average', running_avg_window = None):
    """

    :param x1:
    :param h1: Land ice height vector. Must be resampled to constant spacing
    :param dx:
    :param type:
    :param running_avg_window:
    :return:
    """
    if filter_type == 'running_average': #
        if not running_avg_window: # if no window size was given, run a default 100 m running average smoother
            running_avg_window = 100
        dx = np.median(x1[1:] - x1[:-1])
        filt = np.ones(int(np.round(running_avg_window / dx)))
        h_filt = (1/len(filt)) * np.convolve(filt, h1, mode = 'same')

    # TASK: add bandpass, highpass
    else:
        h_filt = h1
    return h_filt

# -------------------------------------------------------------------------------------------

def differentiate(x_in, h_in):
    """

    :param x_in:
    :param h_in: Land ice height vector. Must be resampled to constant spacing
    :return:
    """
    dh = np.gradient(h_in, x_in)
    return dh

# -------------------------------------------------------------------------------------------
def time_diff(D1,D2):
    r"""Get the time difference between two cycles; units of days

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
    # time step to julian days
    dt = np.abs((t2 - t1).jd)

    return dt

# -------------------------------------------------------------------------------------------

def calculate_velocities(rgt, x_atc, h_li_raw, h_li_diff, lats, lons, segment_ids, times, beams, cycle1, cycle2, product, segment_length, search_width, along_track_step, max_percent_nans, dx, saving = False, write_out_path = '.', prepend = '', spatial_extent = None, map_data_root=None):
    """

    :param x_atc:
    :param h_li_diff:
    :param lats:
    :param lons:
    :param segment_ids:
    :param beams:
    :param cycle1:
    :param cycle2:
    :param segment_length:
    :param search_width:
    :param along_track_step:
    :param max_percent_nans:
    :param dx:
    :return:
    """
    # takes output of load_data_by_rgt

    from scipy.signal import correlate, detrend
    import glob

    ### Create dictionaries to put info in
    velocities = {}
    correlations = {}
    lags = {}
    x_atcs_for_velocities = {}
    latitudes = {}
    longitudes = {}
    midpoints_x_atc = {}
    midpoints_lats = {}
    midpoints_lons = {}
    midpoints_seg_ids = {}
    midpoints_xy = {}

    # inside rgt loop
    velocities[rgt] = {}
    correlations[rgt] = {}
    lags[rgt] = {}
    midpoints_x_atc[rgt] = {}
    midpoints_lats[rgt] = {}
    midpoints_lons[rgt] = {}
    midpoints_seg_ids[rgt] = {}
    midpoints_xy[rgt] = {}

    ### Determing Elapsed time between cycles
    t1_string = times[cycle1]['gt1l'][0].astype(str)  # figure out later if just picking hte first one it ok
    t1 = Time(t1_string)
    t2_string = times[cycle2]['gt1l'][0].astype(str)  # figure out later if just picking hte first one it ok
    t2 = Time(t2_string)
    dt = (t2 - t1).jd  # difference in julian days

    if saving:
        ### Where to save the results:
        h5_file_out = write_out_path + prepend + '_rgt' + rgt + '.hdf5'

        ### Save some metadata
        with h5py.File(h5_file_out, 'w') as f:
            f['dx'] = dx
            f['product'] = product
            f['segment_length'] = segment_length
            f['search_width'] = search_width
            f['along_track_step'] = along_track_step
            f['max_percent_nans'] = max_percent_nans
            # f['smoothing'] = smoothing
            # f['smoothing_window_size'] = smoothing_window_size
            f['process_date'] = str(Time.now().value)
            f['rgt'] = rgt


    ### Loop over each beam
    for beam in beams:

        ### Determine x1s, which are the x_atc coordinates at which each cut out window begins
        # To be common between both repeats, the first point x1 needs to be the larger first value between repeats
        min_x_atc_cycle1 = x_atc[cycle1][beam][0]
        min_x_atc_cycle2 = x_atc[cycle2][beam][0]

        # pick out the track that starts at greater x_atc, and use that as x1s vector
        if min_x_atc_cycle1 != min_x_atc_cycle2:
            x1 = np.nanmax([min_x_atc_cycle1, min_x_atc_cycle2])
            cycle_n = np.arange(0, 2)[[min_x_atc_cycle1, min_x_atc_cycle2] == x1][0]
            if cycle_n == 0:
                cycletmp = cycle2
            elif cycle_n == 1:
                cycletmp = cycle1
            # n_segments_this_track = (len(x_atc[cycletmp][beam]) - 2 * search_width/dx) / (along_track_step/dx)

            ### Generate the x1s vector, in the case that the repeat tracks don't start in the same place
            x1s = x_atc[cycletmp][beam][
                  int(search_width / dx) + 1:-int(segment_length / dx) - int(search_width / dx):int(along_track_step / dx)]
            #### ! this line updated 2020 07 23
            # x1s = x_atc[cycletmp][beam][int(search_width/dx)+1::int(search_width/dx)]
            # start at search_width/dx in, so the code never tries to get data outside the edges of this rgt
            # add 1 bc the data are differentiated, and h_li_diff is therefore one point shorter

        elif min_x_atc_cycle1 == min_x_atc_cycle2:  # doesn't matter which cycle
            ### Generate the x1s vector, in the case that the repeat tracks do start in the same place
            x1s = x_atc[cycle1][beam][
                  int(search_width / dx) + 1:-int(segment_length / dx) - int(search_width / dx):int(along_track_step / dx)]
            #### ! this line updated 2020 07 23
            # x1s = x_atc[cycle1][beam][int(search_width/dx)+1::int(search_width/dx)]

        ### Determine xend, where the x1s vector ends: smaller value for both beams, if different
        max_x_atc_cycle1 = x_atc[cycle1][beam][-1] - search_width / dx
        max_x_atc_cycle2 = x_atc[cycle2][beam][-1] - search_width / dx
        smallest_xatc = np.min([max_x_atc_cycle1, max_x_atc_cycle2])
        ixmax = np.where(x1s >= (smallest_xatc - search_width / dx))
        if len(ixmax[0]) >= 1:
            ixtmp = ixmax[0][0]
            x1s = x1s[:ixtmp]

        ### Create vectors to store results in
        velocities[rgt][beam] = np.empty_like(x1s)
        correlations[rgt][beam] = np.empty_like(x1s)
        lags[rgt][beam] = np.empty_like(x1s)
        midpoints_x_atc[rgt][beam] = np.empty_like(x1s)
        midpoints_lats[rgt][beam] = np.empty_like(x1s)
        midpoints_lons[rgt][beam] = np.empty_like(x1s)
        midpoints_seg_ids[rgt][beam] = np.empty_like(x1s)

        # midpoints_x_atc = np.empty(np.shape(x1s))  # for writing out
        # midpoints_lat = np.empty(np.shape(x1s))  # for writing out
        # midpoints_lon = np.empty(np.shape(x1s))  # for writing out
        # midpoints_seg_ids = np.empty(np.shape(x1s))  # for writing out

        ### Entire x_atc vectors for both cycles
        x_full_t1 = x_atc[cycle1][beam]
        x_full_t2 = x_atc[cycle2][beam]

        ### Loop over x1s, positions along track that each window starts at
        for xi, x1 in enumerate(x1s):

            ### Cut out data: small chunk of data at time t1 (first cycle)
            ix_x1 = np.arange(len(x_full_t1))[x_full_t1 >= x1][0]  # Index of first point that is greater than x1
            ix_x2 = ix_x1 + int(
                np.round(segment_length / dx))  # ix_x1 + number of datapoints within the desired segment length
            x_t1 = x_full_t1[ix_x1:ix_x2]  # cut out x_atc values, first cycle
            lats_t1 = lats[cycle1][beam][ix_x1:ix_x2]  # cut out latitude values, first cycle
            lons_t1 = lons[cycle1][beam][ix_x1:ix_x2]  # cut out longitude values, first cycle
            seg_ids_t1 = segment_ids[cycle1][beam][ix_x1:ix_x2]  # cut out segment_ids, first cycle
            h_li1 = h_li_diff[cycle1][beam][
                    ix_x1 - 1:ix_x2 - 1]  # cut out land ice height values, first cycle; start 1 index earlier because
            # the h_li_diff data are differentiated, and therefore one sample shorter

            # Find segment midpoints; this is the position where we will assign the velocity measurement from each window
            n = len(x_t1)
            midpt_ix = int(np.floor(n / 2))
            midpoints_x_atc[rgt][beam][xi] = x_t1[midpt_ix]
            midpoints_lats[rgt][beam][xi] = lats_t1[midpt_ix]
            midpoints_lons[rgt][beam][xi] = lons_t1[midpt_ix]
            midpoints_seg_ids[rgt][beam][xi] = seg_ids_t1[midpt_ix]

            ### Cut out data: wider chunk of data at time t2 (second cycle)
            ix_x3 = ix_x1 - int(np.round(search_width / dx))  # extend on earlier end by number of indices in search_width
            ix_x4 = ix_x2 + int(np.round(search_width / dx))  # extend on later end by number of indices in search_width
            x_t2 = x_full_t2[ix_x3:ix_x4]  # cut out x_atc values, second cycle
            h_li2 = h_li_diff[cycle2][beam][
                    ix_x3 - 1:ix_x4 - 1]  # cut out land ice height values, second cycle; start 1 index earlier because
            # the h_li_diff data are differentiated, and therefore one sample shorter

            ### Determine number of nans in each data chunk
            n_nans1 = np.sum(np.isnan(h_li_raw[cycle1][beam][ix_x1:ix_x2]))
            n_nans2 = np.sum(np.isnan(h_li_raw[cycle2][beam][ix_x3:ix_x4]))

            ### Only process if there are fewer than 10% nans in either data chunk:
            if (n_nans1 / len(h_li1) <= max_percent_nans / 100) and (n_nans2 / len(h_li2) <= max_percent_nans / 100):

                # Detrend both chunks of data
                h_li1 = detrend(h_li1, type='linear')
                h_li2 = detrend(h_li2, type='linear')

                # Normalize both chunks of data, if desired
                # h_li1 = h_li1 / np.nanmax(np.abs(h_li1))
                # h_li2 = h_li2 / np.nanmax(np.abs(h_li2))

                ### Correlate the old and new data
                # We made the second data vector longer than the first, so the valid method returns values
                corr = correlate(h_li1, h_li2, mode='valid', method='direct')

                ### Normalize correlation function by autocorrelations
                # Normalizing coefficient changes with each step along track; this section determines a changing along track normalizing coefficiant
                coeff_a_val = np.sum(h_li1 ** 2)
                coeff_b_val = np.zeros(len(h_li2) - len(h_li1) + 1)
                for shift in range(len(h_li2) - len(h_li1) + 1):
                    h_li2_section = h_li2[shift:shift + len(h_li1)]
                    coeff_b_val[shift] = np.sum(h_li2_section ** 2)
                norm_vec = np.sqrt(coeff_a_val * coeff_b_val)
                corr_normed = corr / np.flip(
                    norm_vec)  # i don't really understand why this has to flip, but otherwise it yields correlation values above 1...
                # TASK: figure out why the flip is needed

                ### Create a vector of lags for the correlation function
                lagvec = np.arange(- int(np.round(search_width / dx)), int(search_width / dx) + 1, 1)  # for mode = 'valid'

                ### Convert lag to distance
                shift_vec = lagvec * dx

                ### ID peak correlation coefficient
                ix_peak = np.arange(len(corr_normed))[corr_normed == np.nanmax(corr_normed)][0]

                ### Save correlation coefficient, best lag, velocity, etc at the location of peak correlation coefficient
                best_lag = lagvec[ix_peak]
                best_shift = shift_vec[ix_peak]

                #                                 ### ID peak correlation coefficient; 'raw' = no sub-sampling
                #                                 plotting = False
                #                                 best_lag, peak_corr_value = find_correlation_peak(lagvec, shift_vec, corr_normed, max_width, min_width, dx_interp, type, plotting)
                #                                 best_shift = best_lag * dx

                velocities[rgt][beam][xi] = best_shift / (dt / 365)
                correlations[rgt][beam][xi] = corr_normed[ix_peak]
                lags[rgt][beam][xi] = lagvec[ix_peak]


            else:
                ### If there are too many nans, just save a nan
                velocities[rgt][beam][xi] = np.nan
                correlations[rgt][beam][xi] = np.nan
                lags[rgt][beam][xi] = np.nan

        ### Output xy locations of midpoints as well
        midpoints_xy[rgt][beam] = np.array(pyproj.Proj(3031)(midpoints_lons[rgt][beam], midpoints_lats[rgt][beam]))
        #extract measures veloc outside of this loop

        if saving:
            measures_Vx_path = glob.glob( map_data_root + '*_Vx.tif')[0]
            measures_Vy_path = glob.glob( map_data_root + '*_Vy.tif')[0]

            ### Add velocities to hdf5 file for each beam
            with h5py.File(h5_file_out, 'a') as f:
                f[beam + '/x_atc'] = midpoints_x_atc[rgt][beam]  # assign x_atc value of half way along the segment
                f[beam + '/latitudes'] = midpoints_lats[rgt][beam]  # assign x_atc value of half way along the segment
                f[beam + '/longitudes'] = midpoints_lons[rgt][beam]  # assign x_atc value of half way along the segment
                f[beam + '/velocities'] = velocities[rgt][beam]  # assign x_atc value of half way along the segment
                f[beam + '/correlation_coefficients'] = correlations[rgt][beam]  # assign x_atc value of half way along the segment
                f[beam + '/best_lags'] = lags[rgt][beam]  # assign x_atc value of half way along the segment
                f[beam + '/segment_ids'] = midpoints_seg_ids[rgt][beam]
                f[beam + '/first_cycle_time'] = str(Time(times[cycle1][beam][0]))
                f[beam + '/second_cycle_time'] = str(Time(times[cycle2][beam][0]))
                f[beam + '/Measures_v_along'] = get_measures_along_track_velocity(midpoints_xy[rgt][beam][0], midpoints_xy[rgt][beam][1], spatial_extent,
                                                                          measures_Vx_path, measures_Vy_path)

    return velocities, correlations, lags, midpoints_x_atc, midpoints_xy, midpoints_lons, midpoints_lats

def velocity_old(x_in1, x_in2, dh1, dh2, dt, segment_length=2000, search_width=1000, along_track_step=100, dx=20, corr_threshold=.65):
    r"""Calculate along-track velocity by correlating the along-track elevation between repeat acquisitions.
    Works for a single set of repeat beams

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
    output_xs : array # REMOVED THIS ONE, calculate inside the loop
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

    # Generate vector of along window starting positions in the along-track vector
    x1s = x_in[int(search_width/dx)+1:-int(segment_length/dx) - int(search_width/dx):int(along_track_step/dx)]
    # (This is the one that was not quite working correctly in our original notebooks)

    # middle point of each window, where we actually assign the velocity
    # xs_middle =


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

