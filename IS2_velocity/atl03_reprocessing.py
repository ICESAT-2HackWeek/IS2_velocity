#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

def reinterpolate_atl03(x_in,h_in,x_step=5,x_win=10):
    r"""Recreate an ATL06-like product at finer resolution.

    The ICESat-2 ATL06 product is created by linear regression
    on the ATL03 photon cloud every 20 meters; then, the ATL06
    elevation is evaluated at the center-point of the regression.
    This funtion is meant to recrate this processing flow but for any
    arbitrary step size.

    Parameters
    ------
    x_in : array
           along-track photon locations
    h_in : array
           along-track photon height
    x_step : float; optional
             along-track step between points in the desired output array
    x_win : float; optional
            window over which the liner fit is done (window is in both directions)

    Output
    ------
    xs : array
         along-track distance for reinterpolated points
    hs : array
         heights corresponding to the distances in 'xs'
    """

    # Set up the output arrays to be filled
    xs = np.arange(np.nanmin(x_in),np.nanmax(x_in),x_step)
    hs = np.empty_like(xs)
    print('Total Length:',len(xs))

    for i,x in enumerate(xs):
        if i%100 == 0:
            print(i)
        # Find all the points within x_win meters of the current x location
        idxs = np.logical_and(x_in>x-x_win,x_in<x+x_win)
        # Set up a linear regression on the points inside the window
        x_fit = x_in[idxs]
        h_fit = h_in[idxs]
        # index the points that are not nans (these will be used in the regression)
        idx_nonan = ~np.isnan(x_fit) & ~np.isnan(h_fit)
        if ~np.any(idx_nonan):
            # If all the points in the window are nans, output nan
            hs[i] = np.nan
        else:
            fit = np.polyfit(x_fit[idx_nonan],h_fit[idx_nonan],1)
            # Evaluate the regression at the centerpoint to get the height
            hs[i] = np.polyval(fit,x)

    return xs,hs
