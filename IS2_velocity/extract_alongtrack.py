#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
try:
    import pointCollection as pc
except:
    print('Continuing withough pointCollection')
from pyproj import Proj
proj_stere = Proj('epsg:3031')

def get_measures(data, spatial_extent, measures_path):
    r"""Add the Measures surface velocity to the local data dictionary.

    Parameters
    ------
    data : dict
        data dictionary
    spatial_extent : list
        bounding box of the interest area in the format
        (e.g. [-65, -86, -55, -81] == [min_lon, min_lat, max_lon, max_lat])
    measures_path : str
        local path to velocity data where Vx and Vy are located (format Vx.tif and Vy.tif)

    Output
    ------
    data : dict
    """

    #fix with if statement about type of list or array DONE
    if type(spatial_extent) == type([]):
        spatial_extent = np.array([spatial_extent])
    lat=spatial_extent[[1, 3, 3, 1, 1]]
    lon=spatial_extent[[2, 2, 0, 0, 2]]

    # project the coordinates to Antarctic polar stereographic
    meas_xy=np.array(proj_stere(lon, lat))
    # get the bounds of the projected coordinates
    XR=[np.nanmin(meas_xy[0,:]), np.nanmax(meas_xy[0,:])]
    YR=[np.nanmin(meas_xy[1,:]), np.nanmax(meas_xy[1,:])]

    Measures_vx=pc.grid.data().from_geotif(measures_path+'Vx.tif', bounds=[XR, YR])
    Measures_vy=pc.grid.data().from_geotif(measures_path+'Vy.tif', bounds=[XR, YR])

    beam1 = list(data['midpoints_lats'].keys())[0]
    lats = data['midpoints_lats'][beam1]
    lons = data['midpoints_lons'][beam1]
    x,y = proj_stere(lons,lats)

    vx = Measures_vx.interp(x,y)
    vy = Measures_vy.interp(x,y)

    #Solve for angle to rotate Vy to be along track and Vx to be across track
    xL = x[1] - x[0]
    yL = y[1] - y[0]
    theta = np.arctan(yL/xL)

    # Do the rotation
    data['meas_v_along'] = vx*np.cos(theta) + vy*np.cos(np.pi/2.-theta)
    data['meas_v_across'] = vx*np.sin(theta) - vy*np.sin(np.pi/2.-theta)

    return data

# =============================================================================
#     import math
#     xL=abs((is2_dict['x'][0])-(is2_dict['x'][1]))
#     yL=abs((is2_dict['y'][0])-(is2_dict['y'][1]))
# 
#     #decides if is descending or ascending path
#     if is2_dict['x'][0]-is2_dict['x'][1] < 0:
# 
#         theta_rad=math.atan(xL/yL)
#         #theta_deg=theta_rad*180/math.pi
#         is2_dict['v_along']=vy/math.cos(theta_rad)
#         is2_dict['v_across']=vx/math.cos(theta_rad)
# 
#     else:
# 
#         theta_rad=math.atan(xL/yL)
#         #theta_deg=theta_rad*180/math.pi
#         is2_dict['v_along']=vy/math.sin(theta_rad)
#         is2_dict['v_across']=vx/math.sin(theta_rad)
# 
#     return is2_dict
# =============================================================================


