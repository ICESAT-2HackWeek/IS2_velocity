#!/usr/bin/env python
# -*- coding: utf-8 -*-

import glob, os, h5py
import numpy as np
from IS2_velocity.correlation_processing import determine_x1s, correlate_single_x1, time_diff
try:
    import pointCollection as pc
except:
    print('Continuing without pointCollection.')
import matplotlib.pyplot as plt
from pyproj import Proj
proj_stere = Proj('epsg:3031')

def plot_measures_along_track_comparison(rgt, beams, data_path, correlation_threshold, spatial_extent, plot_out_location, map_data_root,
                                         velocity_number, flip_meas_sign=False, close=False):
    """

    :param rgt:
    :param out_path: path to folder where data written out by the correlation step is saved; data is re-loaded in this step
    :param correlation_threshold:
    :param plot_out_location:
    :param map_data_root:
    :param velocity_number:
    :param close:
    :return:
    """ # currently just the first velocity determination, veloc0
    # out_path is where the xcorr results are stored
    # plot_out_location is where to save the plot
    # map_data_root is where the map data are stored, specifically must contain moa_2009_1km.tif for this specific code to work

    file = glob.glob(data_path + '*_rgt' + rgt + '.hdf5')[0]

    # < eventually, make velocity number two digits
    # '_veloc' + str(velocity_number).zfill(2) +

    if glob.glob(file):

        plt.close('all')
        fig = plt.figure(figsize=[11, 8])
        grid = plt.GridSpec(6, 2, wspace=0.4, hspace=0.3)

        try:
            ### MOA parameters
            moa_datapath = map_data_root  # '/srv/tutorial-data/land_ice_applications/'

            lat = spatial_extent[[1, 3, 3, 1, 1]]
            lon = spatial_extent[[2, 2, 0, 0, 2]]

            # project the coordinates to Antarctic polar stereographic
            xy = np.array(proj_stere(lon, lat))
            # get the bounds of the projected coordinates
            XR = [np.nanmin(xy[0, :]), np.nanmax(xy[0, :])]
            YR = [np.nanmin(xy[1, :]), np.nanmax(xy[1, :])]
            #             MOA=pc.grid.data().from_geotif(os.path.join(moa_datapath, 'MOA','moa_2009_1km.tif'), bounds=[XR, YR])

            MOA = pc.grid.data().from_geotif(os.path.join(moa_datapath, 'moa_2009_1km.tif'), bounds=[XR, YR])

            haxMOA = fig.add_subplot(grid[0:4, 1])
            MOA.show(ax=haxMOA, cmap='gray', clim=[14000, 17000])
        except:
            pass

        with h5py.File(file, 'r') as f:
            beam = beams[0]
            lats = f[f'/{beam}/latitudes'][()]
            lons = f[f'/{beam}/longitudes'][()]
            meas_v = f[f'/{beam}/Measures_v_along'][()]
            meas_xy = np.array(proj_stere(lons, lats))
            if flip_meas_sign:
                meas_v *= -1

            for ib, beam in enumerate(beams):
                hax0 = fig.add_subplot(grid[ib, 0])
                # 1hax1=fig.add_subplot(212)
                # hax1.set_title('measures ' )
                if ib == 0:
                    hax0.set_title('velocs vs measures ' + rgt)

                lats = f[f'/{beam}/latitudes'][()]
                lons = f[f'/{beam}/longitudes'][()]
                coeffs = f[f'/{beam}/correlation_coefficients'][()]
                velocs = f[f'/{beam}/velocities'][()]
                xy = np.array(proj_stere(lons, lats))

                ixs0 = coeffs <= correlation_threshold
                ixs = coeffs > correlation_threshold

                h0 = hax0.scatter(xy[0], velocs, 1, coeffs, vmin=0, vmax=1, cmap='viridis')
                h1 = hax0.plot(meas_xy[0], meas_v, 'k-')
                # whether v_along is + or - must depend on ascending vs descending; not done correctly yet

                hax0.set_ylim(-800, 800)
                c = plt.colorbar(h0, ax=hax0)
                c.set_label('Correlation\ncoefficient\n(-1 -> 1)')

                try:
                    h2 = haxMOA.scatter(xy[0][ixs0], xy[1][ixs0], 0.02, 'k')
                    h3 = haxMOA.scatter(xy[0][ixs], xy[1][ixs], 0.15, velocs[ixs], vmin=-800, vmax=800, cmap='plasma')
                    if ib == 0:
                        c = plt.colorbar(h3, ax=haxMOA)
                        c.set_label('Along-track velocity (m/yr)')
                except:
                    print('error')
                    pass

        outfile = plot_out_location + 'rgt' + rgt + '.' + beam + '_vs_measures_veloc' + str(velocity_number).zfill(
            2) + '.png'
        plt.savefig(outfile, dpi=200)
        if close == True:
            plt.close('all')

def plot_correlation_function(x1, data, cycle1, cycle2, beam, search_width=1000, segment_length=2000, along_track_step=100,
                                    max_percent_nans=10, dx=20):
    r""" Re-calculate correlation function for a particular along-track starting place.
    Useful for visualizing quality of correlation functions

    Parameters
    ------
    :param data:
        Dictionary containing loaded hdf5 data and correlation processing results
    :param cycle1:
        Name of earlier cycle to compare
    :param cycle2:
        Name of later cycle to compare
    :param beam:
        Name of specific beam to compare
    :param xi:
        Along-track window number

    Output
    ------
    fig
    """

    ### This function reproduces the correlation processing, but only for a single value of along-track distance

    ### Determine x1s, the along-track distances at which to start each window
    x1s = determine_x1s(data, cycle1, cycle2, beam, search_width=search_width, segment_length=segment_length,
                                    along_track_step=along_track_step, dx=dx)

    ### Determine xi, the index of the current along track distance in the x1s vector
    xi = np.where(x1s == x1)[0][0]

    ### Correlate for that single value of x1
    data_tmp = correlate_single_x1(x1, data, cycle1, cycle2, beam, search_width=search_width, segment_length=segment_length,
                               max_percent_nans=max_percent_nans, dx=dx, return_correlation_function=True)


    ### Cut out data; code taken from correlate_single_x1s
    npts_short = int( np.round(segment_length / dx)) # how many data points in the later, shorter vector to correlate
    npts_extra_long = int(np.round(search_width / dx)) # how many extra datapoints before/after in the earlier, longer vector to correlate

    # x_full_t1 = data['x_atc'][cycle1][beam]

    ### Cut out data, t1 (first cycle)
    ix_x1 = np.arange(len(data['x_atc'][cycle1][beam]))[data['x_atc'][cycle1][beam] >= x1][
        0]  # Index of first point that is greater than x1
    ix_x2 = ix_x1 + npts_short  # ix_x1 + number of datapoints within the desired segment length
    # x_t1 = x_full_t1[ix_x1:ix_x2]  # cut out x_atc values, first cycle
    # lats_t1 = data['latitude'][cycle1][beam][ix_x1:ix_x2]  # cut out latitude values, first cycle
    # lons_t1 = data['longitude'][cycle1][beam][ix_x1:ix_x2]  # cut out longitude values, first cycle
    # seg_ids_t1 = data['segment_ids'][cycle1][beam][ix_x1:ix_x2]  # cut out segment_ids, first cycle
    x_atc1 = data['x_atc'][cycle1][beam][ix_x1:ix_x2]
    h_li1 = data['h_li'][cycle1][beam][ix_x1:ix_x2]
    h_li_diff1 = data['h_li_diff'][cycle1][beam][ ix_x1 - 1:ix_x2 - 1]  # cut out land ice height values, first cycle; start 1 index earlier because
    # the h_li_diff data are differentiated, and therefore one sample shorter

    ### Cut out data: wider chunk of data at time t2 (second cycle)
    ix_x3 = ix_x1 - npts_extra_long  # extend on earlier end by number of indices in search_width
    ix_x4 = ix_x2 + npts_extra_long  # extend on later end by number of indices in search_width
    x_atc2 = data['x_atc'][cycle1][beam][ix_x3:ix_x4]
    h_li2 = data['h_li'][cycle2][beam][ix_x3:ix_x4]  # cut out land ice height values, second cycle; start 1 index earlier because
    h_li_diff2 = data['h_li_diff'][cycle2][beam][ ix_x3 - 1:ix_x4 - 1]  # cut out land ice height values, second cycle; start 1 index earlier because
    # the h_li_diff data are differentiated, and therefore one sample shorter

    veloc = data_tmp['velocities'][beam][xi]
    correlation_value = data_tmp['correlations'][beam][xi]
    best_lag = data_tmp['lags'][beam][xi]

    shift_vec = data['correlation_functions'][beam][xi]['shift_vec'] # units of meters
    correlation_function = data['correlation_functions'][beam][xi]['correlation_function']

      ### Create figure
    fig = plt.figure(figsize=[11, 8])
    grid = plt.GridSpec(5, 1, wspace=0.4, hspace=0.3)

    ### Plot raw data for this comparison
    ax0 = fig.add_subplot(grid[0, 0])
    cycle2_plot, = ax0.plot(x_atc2 - x_atc1[0], h_li2, 'r')
    cycle1_plot, = ax0.plot(x_atc1 - x_atc1[0], h_li1, 'k')
    ax0.legend((cycle1_plot, cycle2_plot), (cycle1, cycle2), loc='upper right')
    ax0.set_ylabel('Ice Surface\nElevation (m)')
    ax0.set_xlim(np.min(x_atc2 - x_atc1[0]), np.max(x_atc2 - x_atc1[0]))
    ax0.text(x_atc2[3] - x_atc1[0], np.max(h_li2)-0.25, 'Raw data')

    ### Plot differentiated data for this comparison
    ax1 = fig.add_subplot(grid[1, 0])
    ax1.plot(x_atc2 - x_atc1[0], h_li_diff2, 'r')
    ax1.plot(x_atc1 - x_atc1[0], h_li_diff1, 'k')
    ax1.set_ylabel('Differentiated\nElevation')
    ax1.text(x_atc2[3] - x_atc1[0], 0.8*np.max(h_li_diff2), 'Raw data, differentiated')
    ax1.set_xlim(np.min(x_atc2 - x_atc1[0]), np.max(x_atc2 - x_atc1[0]))


    ### Plot shifted data for this comparison
    ax2 = fig.add_subplot(grid[2, 0])
    ax2.plot(x_atc2 - x_atc1[0], h_li_diff2, 'r')
    ax2.plot(x_atc1 - x_atc1[0] - best_lag * dx, h_li_diff1,
             'k')  # best_lag * dx in case best_lag is sub-sample and can't offset by index
    ax2.set_xlabel('Meters along track')
    ax2.set_ylabel('Differentiated\nElevation')
    ax2.text(x_atc2[3] - x_atc1[0], 0.8*np.max(h_li_diff2), 'Shifted by best lag')
    ax2.set_xlim(np.min(x_atc2 - x_atc1[0]), np.max(x_atc2 - x_atc1[0]))

    ### Plot correlation function for this comparison
    ax3 = fig.add_subplot(grid[4, 0])
    corr_vec_plot, = ax3.plot(shift_vec, correlation_function, 'k')
    ax3.set_xlabel('Shift (meters)')
    ax3.set_ylabel('Correlation coefficient\n(-1 to 1)')
    ax3.set_xlim(np.min(shift_vec), np.max(shift_vec))

    ax3.legend(([corr_vec_plot]), (['Correlation\nFunction']), loc='upper right')

    lag_ix = int(search_width / dx + best_lag)
    ax3.plot(shift_vec[lag_ix], correlation_value, 'ro', markerfacecolor='none')

    dt = time_diff(data_tmp, cycle1, cycle2, beam)

    info = 'Peak correlation value of ' + str(np.round(correlation_value,3)) + '\n' \
        'Peak correlation value at shift of ' + str(shift_vec[lag_ix]) +\
            ' m\nCorresponding Velocity = ' + str(np.round(veloc,1)) + ' m/yr'
    ax3.set_title(info)

    info2 = 'x1 = ' + str(np.round(x1,3)) + '; xi = ' + str(xi) + \
        '\nCycles: ' + ','.join([cycle1, cycle2]) + \
        '\nBeam: ' + beam
    fig.suptitle(info2)

    return fig

