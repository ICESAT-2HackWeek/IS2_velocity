#!/usr/bin/env python
# -*- coding: utf-8 -*-

from astropy.time import Time
import h5py

def save_is2_velocity(data,rgt,write_out_path,prepend,product,map_data_root,
                    cycle1,cycle2,dx,segment_length,search_width,along_track_step,
                    max_percent_nans,spatial_extent):
    """Save the dictionary as an hdf5 file.

    Parameters
    ------
    data : dict
           data dictionary
    rgt
    write_out_path
    prepend
    product
    map_data_root
    cycle1
    cycle2
    dx
    segment_length
    search_width
    along_track_step
    max_percent_nans
    spatial_extent

    Output
    ------
    """

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

    for beam in data['velocities'].keys():

        ### Add velocities to hdf5 file for each beam
        with h5py.File(h5_file_out, 'a') as f:
            f[beam + '/x_atc'] = data['midpoints_x_atc'][beam]  # assign x_atc value of half way along the segment
            f[beam + '/latitudes'] = data['midpoints_lats'][beam]  # assign x_atc value of half way along the segment
            f[beam + '/longitudes'] = data['midpoints_lons'][beam]  # assign x_atc value of half way along the segment
            f[beam + '/velocities'] = data['velocities'][beam]  # assign x_atc value of half way along the segment
            f[beam + '/correlation_coefficients'] = data['correlations'][beam]  # assign x_atc value of half way along the segment
            f[beam + '/best_lags'] = data['lags'][beam]  # assign x_atc value of half way along the segment
            f[beam + '/segment_ids'] = data['midpoints_seg_ids'][beam]
            f[beam + '/first_cycle_time'] = str(Time(data['times'][cycle1][beam][0]))
            f[beam + '/second_cycle_time'] = str(Time(data['times'][cycle2][beam][0]))
            f[beam + '/Measures_v_along'] = data['meas_v_along']
            f[beam + '/Measures_v_across'] = data['meas_v_across']
