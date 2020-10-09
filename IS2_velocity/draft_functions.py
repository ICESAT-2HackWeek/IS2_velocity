#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

### UNFINISHED FUNCTION - basic structure is good
### working on debugging problems with peaks at edges for the sub-sample routine

def find_correlation_peak(lagvec, shift_vec, corr_normed, max_width, min_width, dx_interp, fit_type, plotting):
    # max_width = how many SAMPLES wide the half peak can be, max, samples
    # min_width = the narrowest the half peak can be, samples
    # dx_interp = interpolated sample spacing, in units of SAMPLES

    from scipy.interpolate import interp1d

    # peak of un-interpolated data
    ix_peak = np.arange(len(corr_normed))[corr_normed == np.nanmax(corr_normed)][0]
    ix_peak = np.where(corr_normed == np.nanmax(corr_normed))[0][0]

    if fit_type == 'raw':
        best_lag = lagvec[ix_peak]
        best_shift = shift_vec[ix_peak]
        peak_corr_value = corr_normed[ix_peak]

    else:
        ### Assume we're doing a sub-sample fit
        # find troughs before and after peak

        corr_diff = corr_normed[1:] - corr_normed[:-1] # to find troughs, zero-crossings

        # cases:
        # ix_peak is too close to the start (ix_peak < min_width)
            # calculate right width by finding trough or max/min, use ix_peak as left width
        if ix_peak < min_width:
            left_width = ix_peak # go to left edge
            # determine right width
            if ix_peak == 0:
                ix_next_trough = np.where(corr_diff[ix_peak:] > 0)[0][0] + ix_peak
            else:
                ix_next_trough = np.where(corr_diff[ix_peak-1:] > 0)[0][0] + ix_peak
            width_next_trough = -(ix_peak - ix_next_trough)
            if width_next_trough < min_width:
                right_width = min_width
            elif width_next_trough > max_width:
                right_width = max_width
            else:
                right_width = width_next_trough

        # ix_peak is too close to the end (len(corr_normed) - ix_peak < min_width)
            # calculate left width by finding trough or max/min, use len(corr_normed) - ix_peak as right width
        elif len(corr_normed) - ix_peak < min_width:
            right_width = len(corr_normed) - ix_peak -1
            #determine left width
            ix_prev_trough = ix_peak - np.where(np.flip(corr_diff[:ix_peak-1]) < 0)[0][0]
            width_prev_trough = ix_peak - ix_prev_trough
            if width_prev_trough < min_width:
                left_width = min_width
            elif width_prev_trough > max_width:
                left_width = max_width
            else:
                left_width = width_prev_trough

        # other
        else:
            # calculate both left and right width
            ix_next_trough = np.where(corr_diff[ix_peak-1:] > 0)[0][0] + ix_peak
            ix_prev_trough = ix_peak - np.where(np.flip(corr_diff[:ix_peak-1]) < 0)[0][0]
            # if width is greater than min, pick smaller width, so is same on both sides
            width = np.min([ix_peak - ix_prev_trough, -(ix_peak - ix_next_trough)])
            if width > max_width:
                width = max_width
            elif width < min_width:
                width = min_width
            left_width = width
            right_width = width

#         # deal with edges; just go to edge
#         if ix_peak + width >= len(corr_normed):
#             right_width = len(corr_normed) - ix_peak -1
#         else:
#             right_width = width

#         if width>= ix_peak:
#             left_width = ix_peak
#         else:
#             left_width = width

        # cut out data and an x vector before and after the peak (for plotting)
        dcut = corr_normed[ix_peak - left_width: ix_peak + right_width +1]
        xvec_cut = np.arange(ix_peak - left_width,  ix_peak + right_width + 1)
        lagvec_cut = lagvec[xvec_cut]

        if fit_type == 'parabola':
            # fit parabola; in these data, usually a poor fit
            fit_coeffs = np.polyfit(xvec_cut, dcut, deg=2, full = True)

            # create xvec to interpolate the parabola onto and compute
            interp_xvec = np.arange(ix_peak - left_width,  ix_peak + right_width + 1, 0.01)
            fit = np.polyval(p = fit_coeffs[0], x = interp_xvec)

            # interp lagvec
            lagvec_interp = np.arange(lagvec[ix_peak - left_width], lagvec[ix_peak + right_width + 1], dx_interp)

            # compute peak
            ix_peak_interp = np.where(fit == np.max(fit))[0][0]
            best_lag = interp_xvec[ix_peak_interp]
            peak_corr_value = fit[ix_peak_interp]

        elif fit_type == 'cubic':
#             # find troughs before and after peak
#             corr_diff = corr_normed[1:] - corr_normed[:-1]
#             ix_next_trough = np.where(corr_diff[ix_peak:] > 0)[0][0] + ix_peak
#             ix_prev_trough = ix_peak - np.where(np.flip(corr_diff[:ix_peak]) < 0)[0][0]
#             width = np.min([ix_peak - ix_prev_trough, -(ix_peak - ix_next_trough)])
#             if width > max_width:
#                 width = max_width
#             elif width < min_width:
#                 width = min_width

#             # cut out data and an x vector before and after the peak (for plotting)
#             dcut = corr_normed[ix_peak - width: ix_peak + width +1]
#             xvec_cut = np.arange(ix_peak - width,  ix_peak + width + 1)
#             lagvec_cut = lagvec[xvec_cut]

#             if np.sum(np.isnan(dcut)) >0:
#                 print(str(np.sum(np.isnan(dcut))) + ' nans')

            # cubic spline interpolation
            # create xvector to interpolate onto
            interp_xvec = np.arange(ix_peak - left_width,  ix_peak + right_width , dx_interp)
            # create interpolation function
            f = interp1d(xvec_cut, dcut, kind = 'cubic')
            # evaluate interpolation function on new xvector
            fit = f(interp_xvec)

            # interpolate lagvec
            lagvec_interp = np.arange(lagvec[ix_peak - left_width], lagvec[ix_peak + right_width], dx_interp)

            # compute peak
            ix_peak_interp = np.where(fit == np.max(fit))[0][0]
            best_lag = lagvec_interp[ix_peak_interp]
            peak_corr_value = fit[ix_peak_interp]

    return best_lag, peak_corr_value
