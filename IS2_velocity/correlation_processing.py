#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy.signal import correlate,detrend
from astropy.time import Time


# -------------------------------------------------------------------------------------------
def time_diff(data,cycle1,cycle2,beam='gt1l'):
    r"""Get the time difference between two cycles; units of days

    Parameters
    ------
    data : dict
        data dictionary
    cycle1 : int
        first cycle to be compared
    cycle2 : int
        second cycle to be compared
    beam : str
        which beam to compare

    Output
    ------
    dt : float
         time step
    """

    # get time strings
    t1_string = data['times'][cycle1][beam]
    t2_string = data['times'][cycle2][beam]
    # use astropy to get time from strings
    t1 = Time(t1_string)
    t2 = Time(t2_string)
    # time step to julian days
    dt = np.abs((t2 - t1).jd)

    return dt

# -------------------------------------------------------------------------------------------

def calculate_velocities(data, cycle1, cycle2, beams, *args, **kwargs):
    r"""Calculate velocities by correlating two cycles.
    Can be done over a list of beams.

    Parameters
    ------
    data : dict
        data dictionary
    cycle1 : int
        first cycle to be compared
    cycle2 : int
        second cycle to be compared
    beams : list
        a list of which beam to compare

    Output
    ------
    data : dict
    """

    ### Create dictionaries to put info in
    variables = ['velocities','correlations','lags','midpoints_x_atc','midpoints_lats','midpoints_lons',
            'midpoints_seg_ids']
    for var in variables:
        data[var] = {}

    ### Determing Elapsed time between cycles
    dt = time_diff(data,cycle1,cycle2)

    ### Loop over each beam
    for beam in beams:
        data  = calculate_velocity_single_beam(data,cycle1,cycle2,beam,dt,*args,**kwargs)

    return data

# -------------------------------------------------------------------------------------------

def calculate_velocity_single_beam(data,cycle1,cycle2,beam,dt,search_width=1000,segment_length=2000,
                                    max_percent_nans=10,along_track_step=100, dx=20):
    r"""Calculate velocities by correlating two cycles.
    Can be done over a list of beams.

    Parameters
    ------
    data : dict
        data dictionary
    cycle1 : int
        first cycle to be compared
    cycle2 : int
        second cycle to be compared
    beams : list
        a list of which beam to compare

    Output
    ------
    data : dict
    """

    ### Determine x1s, which are the x_atc coordinates at which each cut out window begins
    # To be common between both repeats, the first point x1 needs to be the larger first value between repeats
    min_x_atc_cycle1 = data['x_atc_full'][cycle1][beam][0]
    min_x_atc_cycle2 = data['x_atc_full'][cycle2][beam][0]

    # pick out the track that starts at greater x_atc, and use that as x1s vector
    if min_x_atc_cycle1 != min_x_atc_cycle2:
        x1 = np.nanmax([min_x_atc_cycle1, min_x_atc_cycle2])
        cycle_n = np.arange(0, 2)[[min_x_atc_cycle1, min_x_atc_cycle2] == x1][0]
        if cycle_n == 0:
            cycletmp = cycle2
        elif cycle_n == 1:
            cycletmp = cycle1

        ### Generate the x1s vector, in the case that the repeat tracks don't start in the same place
        x1s = data['x_atc_full'][cycletmp][beam][
              int(search_width / dx) + 1:-int(segment_length / dx) - int(search_width / dx):int(along_track_step / dx)]

    elif min_x_atc_cycle1 == min_x_atc_cycle2:  # doesn't matter which cycle
        ### Generate the x1s vector, in the case that the repeat tracks do start in the same place
        x1s = data['x_atc_full'][cycle1][beam][
              int(search_width / dx) + 1:-int(segment_length / dx) - int(search_width / dx):int(along_track_step / dx)]

    ### Determine xend, where the x1s vector ends: smaller value for both beams, if different
    max_x_atc_cycle1 = data['x_atc_full'][cycle1][beam][-1] - search_width / dx
    max_x_atc_cycle2 = data['x_atc_full'][cycle2][beam][-1] - search_width / dx
    smallest_xatc = np.min([max_x_atc_cycle1, max_x_atc_cycle2])
    ixmax = np.where(x1s >= (smallest_xatc - search_width / dx))
    if len(ixmax[0]) >= 1:
        ixtmp = ixmax[0][0]
        x1s = x1s[:ixtmp]

    ### Create vectors to store results in
    variables = ['velocities','correlations','lags','midpoints_x_atc','midpoints_lats','midpoints_lons',
            'midpoints_seg_ids']
    for var in variables:
        data[var][beam] = np.empty_like(x1s)

    ### Entire x_atc vectors for both cycles
    x_full_t1 = data['x_atc_full'][cycle1][beam]

    ### Loop over x1s, positions along track that each window starts at
    for xi, x1 in enumerate(x1s):
        ### Cut out data: small chunk of data at time t1 (first cycle)
        ix_x1 = np.arange(len(data['x_atc_full'][cycle1][beam]))[data['x_atc_full'][cycle1][beam] >= x1][0]  # Index of first point that is greater than x1
        ix_x2 = ix_x1 + int(
            np.round(segment_length / dx))  # ix_x1 + number of datapoints within the desired segment length
        x_t1 = x_full_t1[ix_x1:ix_x2]  # cut out x_atc values, first cycle
        lats_t1 = data['latitude'][cycle1][beam][ix_x1:ix_x2]  # cut out latitude values, first cycle
        lons_t1 = data['longitude'][cycle1][beam][ix_x1:ix_x2]  # cut out longitude values, first cycle
        seg_ids_t1 = data['segment_ids'][cycle1][beam][ix_x1:ix_x2]  # cut out segment_ids, first cycle
        h_li1 = data['h_li_diff'][cycle1][beam][
                ix_x1 - 1:ix_x2 - 1]  # cut out land ice height values, first cycle; start 1 index earlier because
        # the h_li_diff data are differentiated, and therefore one sample shorter

        ### Cut out data: wider chunk of data at time t2 (second cycle)
        ix_x3 = ix_x1 - int(np.round(search_width / dx))  # extend on earlier end by number of indices in search_width
        ix_x4 = ix_x2 + int(np.round(search_width / dx))  # extend on later end by number of indices in search_width
        h_li2 = data['h_li_diff'][cycle2][beam][
                ix_x3 - 1:ix_x4 - 1]  # cut out land ice height values, second cycle; start 1 index earlier because
        # the h_li_diff data are differentiated, and therefore one sample shorter

        # Find segment midpoints; this is the position where we will assign the velocity measurement from each window
        n = len(x_t1)
        midpt_ix = int(np.floor(n / 2))
        data['midpoints_x_atc'][beam][xi] = x_t1[midpt_ix]
        data['midpoints_lats'][beam][xi] = lats_t1[midpt_ix]
        data['midpoints_lons'][beam][xi] = lons_t1[midpt_ix]
        data['midpoints_seg_ids'][beam][xi] = seg_ids_t1[midpt_ix]

        ### Determine number of nans in each data chunk
        n_nans1 = np.sum(np.isnan(data['h_li'][cycle1][beam][ix_x1:ix_x2]))
        n_nans2 = np.sum(np.isnan(data['h_li'][cycle2][beam][ix_x3:ix_x4]))

        ### Only process if there are fewer than 10% nans in either data chunk:
        if not (n_nans1 / len(h_li1) <= max_percent_nans / 100) and (n_nans2 / len(h_li2) <= max_percent_nans / 100):
            ### If there are too many nans, just save a nan
            data['velocities'][beam][xi] = np.nan
            data['correlations'][beam][xi] = np.nan
            data['lags'][beam][xi] = np.nan
        else:
            # Detrend both chunks of data
            h_li1 = detrend(h_li1, type='linear')
            h_li2 = detrend(h_li2, type='linear')

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
            data['velocities'][beam][xi] = shift_vec[ix_peak] / (dt / 365)
            data['correlations'][beam][xi] = corr_normed[ix_peak]
            data['lags'][beam][xi] = lagvec[ix_peak]

    return data
