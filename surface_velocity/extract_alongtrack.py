#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pointCollection as pc
import os, pyproj

#From Ben Smith's code loading in .tif file, running into issues likely with directories
#data_root='/srv/shared/surface_velocity/FIS_Velocity/'
#spatial_extent = np.array([-102, -76, -98, -74.5])
spatial_extent = np.array([-65, -86, -55, -81])
lat=spatial_extent[[1, 3, 3, 1, 1]]
lon=spatial_extent[[2, 2, 0, 0, 2]]
print(lat)
print(lon)
# project the coordinates to Antarctic polar stereographic
xy=np.array(pyproj.Proj(3031)(lon, lat))
# get the bounds of the projected coordinates
XR=[np.nanmin(xy[0,:]), np.nanmax(xy[0,:])]
YR=[np.nanmin(xy[1,:]), np.nanmax(xy[1,:])]
#Originally tried to load data from a local directory, should change to shared directory
#Measures_vx=pc.grid.data().from_geotif(os.path.join(data_root,'Measures2_FIS_Vx.tif'), bounds=[XR, YR])
#Measures_vy=pc.grid.data().from_geotif(os.path.join(data_root,'Measures2_FIS_Vy.tif'), bounds=[XR, YR])

#Rodrigo Gomez Fell computer path @ UC
data_root = '/mnt/user1/Antarctica/Quantarctica3/Glaciology/MEaSUREs Ice Flow Velocity/'

Measures_vx=pc.grid.data().from_geotif(os.path.join(data_root,'anta_phase_map_VX.tif'), bounds=[XR, YR])
Measures_vy=pc.grid.data().from_geotif(os.path.join(data_root,'anta_phase_map_VY.tif'), bounds=[XR, YR])


def add_surface_velocity_to_is2_dict(is2_dict, spatial_extent, path, vel_x, vel_y ):
    """

    is2_dict: Python dictionary with ATL06 track data
    spatial_extent: bounding box of the interest area in the format:
                    (e.g. [-65, -86, -55, -81] == [min_lon, min_lat, max_lon, max_lat])
    path: local path to velocity data
    vel_x: tif velocity raster with x component
    vel_y: tif velocity raster with y component

    """
    data_root = path

    spatial_extent = np.array([spatial_extent])
    lat=spatial_extent[[1, 3, 3, 1, 1]]
    lon=spatial_extent[[2, 2, 0, 0, 2]]
    print(lat)
    print(lon)
    # project the coordinates to Antarctic polar stereographic
    xy=np.array(pyproj.Proj(3031)(lon, lat))
    # get the bounds of the projected coordinates
    XR=[np.nanmin(xy[0,:]), np.nanmax(xy[0,:])]
    YR=[np.nanmin(xy[1,:]), np.nanmax(xy[1,:])]

    Measures_vx=pc.grid.data().from_geotif(os.path.join(data_root,vel_x), bounds=[XR, YR])
    Measures_vy=pc.grid.data().from_geotif(os.path.join(data_root,vel_y), bounds=[XR, YR])

    vx = Measures_vx.interp(is2_dict['x'],is2_dict['y'])
    vy = Measures_vy.interp(is2_dict['x'],is2_dict['y'])

    #Solve for angle to rotate Vy to be along track and Vx to be across track
    import math
    xL=abs((is2_dict['x'][0])-(is2_dict['x'][1]))
    yL=abs((is2_dict['y'][0])-(is2_dict['y'][1]))

    #decides if is descending or ascending path
    if is2_dict['x'][0]-is2_dict['x'][1] < 0:

        theta_rad=math.atan(xL/yL)
        #theta_deg=theta_rad*180/math.pi
        is2_dict['v_along']=vy/math.cos(theta_rad)
        is2_dict['v_across']=vx/math.cos(theta_rad)

    else:

        theta_rad=math.atan(xL/yL)
        #theta_deg=theta_rad*180/math.pi
        is2_dict['v_along']=vy/math.sin(theta_rad)
        is2_dict['v_across']=vx/math.sin(theta_rad)

    return is2_dict


def get_measures_along_track_velocity(x_ps_beam, y_ps_beam , spatial_extent, vel_x_path, vel_y_path):
    """

    is2_dict: Python dictionary with ATL06 track data
    spatial_extent: bounding box of the interest area in the format:
                    (e.g. [-65, -86, -55, -81] == [min_lon, min_lat, max_lon, max_lat])
    path: local path to velocity data
    vel_x: tif velocity raster with x component
    vel_y: tif velocity raster with y component

    """

    #fix with if statement about type of list or array DONE

    if type(spatial_extent) == type([]):

        spatial_extent = np.array([spatial_extent])


    lat=spatial_extent[[1, 3, 3, 1, 1]]
    lon=spatial_extent[[2, 2, 0, 0, 2]]

    # project the coordinates to Antarctic polar stereographic
    xy=np.array(pyproj.Proj(3031)(lon, lat))
    # get the bounds of the projected coordinates
    XR=[np.nanmin(xy[0,:]), np.nanmax(xy[0,:])]
    YR=[np.nanmin(xy[1,:]), np.nanmax(xy[1,:])]

    #Measures_vx=pc.grid.data().from_geotif(os.path.join(data_root,vel_x), bounds=[XR, YR])
    #Measures_vy=pc.grid.data().from_geotif(os.path.join(data_root,vel_y), bounds=[XR, YR])

    Measures_vx=pc.grid.data().from_geotif(vel_x_path, bounds=[XR, YR])
    Measures_vy=pc.grid.data().from_geotif(vel_y_path, bounds=[XR, YR])

    vx = Measures_vx.interp(x_ps_beam,y_ps_beam)
    vy = Measures_vy.interp(x_ps_beam,y_ps_beam)

    #Solve for angle to rotate Vy to be along track and Vx to be across track
    import math
    xL=abs((x_ps_beam[0])-(x_ps_beam[1]))
    yL=abs((y_ps_beam[0])-(y_ps_beam[1]))

    #decides if is descending or ascending path
    if x_ps_beam[0]-x_ps_beam[1] < 0:

        theta_rad=math.atan(xL/yL)
        #theta_deg=theta_rad*180/math.pi
        v_along=vy/math.cos(theta_rad)
        #v_across=vx/math.cos(theta_rad)

    else:

        theta_rad=math.atan(xL/yL)
        #theta_deg=theta_rad*180/math.pi
        v_along=vy/math.sin(theta_rad)
        #v_across=vx/math.sin(theta_rad)

    #Vdiff=vy-v_along
    return v_along

