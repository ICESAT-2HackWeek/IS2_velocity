#!/usr/bin/env python
# -*- coding: utf-8 -*-

import glob, pyproj, os, h5py
import numpy as np
try:
    import pointCollection as pc
except:
    print('Continuing without pointCollection.')
import matplotlib.pyplot as plt

def plot_measures_along_track_comparison(rgt, beams, out_path, correlation_threshold, spatial_extent, plot_out_location, map_data_root,
                                         velocity_number, close=False):
    """

    :param rgt:
    :param out_path: path to folder where data written out by the correlation step is saved; data is re-loaded in this step
    :param correlation_threshold:
    :param plot_out_location:
    :param map_data_root:
    :param velocity_number:
    :param close:
    :return:
    """
    # currently just the first velocity determination, veloc0
    # out_path is where the xcorr results are stored
    # plot_out_location is where to save the plot
    # map_data_root is where the map data are stored, specifically must contain moa_2009_1km.tif for this specific code to work

    file = out_path + 'rgt' + rgt + '_veloc' + str(velocity_number).zfill(
        2) + '.hdf5'  # < eventually, make velocity number two digits

    if glob.glob(file):

        ### MOA parameters
        moa_datapath = map_data_root  # '/srv/tutorial-data/land_ice_applications/'
        # spatial_extent = np.array([-102, -76, -98, -74.5])
        # spatial_extent = np.array([-65, -86, -55, -81])

        lat = spatial_extent[[1, 3, 3, 1, 1]]
        lon = spatial_extent[[2, 2, 0, 0, 2]]
        # project the coordinates to Antarctic polar stereographic
        xy = np.array(pyproj.Proj(3031)(lon, lat))
        # get the bounds of the projected coordinates
        XR = [np.nanmin(xy[0, :]), np.nanmax(xy[0, :])]
        YR = [np.nanmin(xy[1, :]), np.nanmax(xy[1, :])]
        MOA = pc.grid.data().from_geotif(os.path.join(moa_datapath, 'moa_2009_1km.tif'), bounds=[XR, YR])
        #             MOA=pc.grid.data().from_geotif(os.path.join(moa_datapath, 'MOA','moa_2009_1km.tif'), bounds=[XR, YR])

        epsg = 3031  # PS?

        plt.close('all')
        fig = plt.figure(figsize=[11, 8])
        grid = plt.GridSpec(6, 2, wspace=0.4, hspace=0.3)
        haxMOA = fig.add_subplot(grid[0:4, 1])
        MOA.show(ax=haxMOA, cmap='gray', clim=[14000, 17000])

        with h5py.File(file, 'r') as f:
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
                v_along = f[f'/{beam}/Measures_v_along'][()]
                xy = np.array(pyproj.proj.Proj(3031)(lons, lats))

                ixs0 = coeffs <= correlation_threshold
                ixs = coeffs > correlation_threshold

                h0 = hax0.scatter(xy[0], velocs, 1, coeffs, vmin=0, vmax=1, cmap='viridis')
                h1 = hax0.plot(xy[0], v_along, 'k-')
                # whether v_along is + or - must depend on ascending vs descending; not done correctly yet

                hax0.set_ylim(-800, 800)
                c = plt.colorbar(h0, ax=hax0)
                c.set_label('Correlation coefficient (0 -> 1)')

                h2 = haxMOA.scatter(xy[0][ixs0], xy[1][ixs0], 0.02, 'k')
                h3 = haxMOA.scatter(xy[0][ixs], xy[1][ixs], 0.15, velocs[ixs], vmin=-800, vmax=800, cmap='plasma')

        c = plt.colorbar(h3, ax=haxMOA)
        c.set_label('Along-track velocity (m/yr)')

        outfile = plot_out_location + 'rgt' + rgt + '.' + beam + '_vs_measures_veloc' + str(velocity_number).zfill(
            2) + '.png'
        plt.savefig(outfile, dpi=200)
        if close == True:
            plt.close('all')
