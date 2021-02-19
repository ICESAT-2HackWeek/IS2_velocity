#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import os,re,h5py

# -----------------------------------------------------------------------------------------------

def read_HDF5_ATL03(FILENAME, ATTRIBUTES=True, VERBOSE=False):
    """
    ### Adapted from a notebook by Tyler Sutterly 6/14/2910 ###

    #-- PURPOSE: read ICESat-2 ATL03 HDF5 data files

    """
    #-- Open the HDF5 file for reading
    fileID = h5py.File(os.path.expanduser(FILENAME), 'r')

    #-- Output HDF5 file information
    if VERBOSE:
        print(fileID.filename)
        print(list(fileID.keys()))

    #-- allocate python dictionaries for ICESat-2 ATL03 variables and attributes
    IS2_atl03_mds = {}
    IS2_atl03_attrs = {} if ATTRIBUTES else None

    #-- read each input beam within the file
    IS2_atl03_beams = [k for k in fileID.keys() if bool(re.match('gt\d[lr]',k))]
    for gtx in IS2_atl03_beams:
        IS2_atl03_mds[gtx] = {}
        IS2_atl03_mds[gtx]['heights'] = {}
        IS2_atl03_mds[gtx]['geolocation'] = {}
    #         IS2_atl03_mds[gtx]['bckgrd_atlas'] = {}
        IS2_atl03_mds[gtx]['geophys_corr'] = {}
        #-- get each HDF5 variable
        #-- ICESat-2 Measurement Group
        for key,val in fileID[gtx]['heights'].items():
            IS2_atl03_mds[gtx]['heights'][key] = val[:]
        #-- ICESat-2 Geolocation Group
        for key,val in fileID[gtx]['geolocation'].items():
            IS2_atl03_mds[gtx]['geolocation'][key] = val[:]
    #         #-- ICESat-2 Background Photon Rate Group
    #         for key,val in fileID[gtx]['bckgrd_atlas'].items():
    #             IS2_atl03_mds[gtx]['bckgrd_atlas'][key] = val[:]
        #-- ICESat-2 Geophysical Corrections Group: Values for tides (ocean,
        #-- solid earth, pole, load, and equilibrium), inverted barometer (IB)
        #-- effects, and range corrections for tropospheric delays
        for key,val in fileID[gtx]['geophys_corr'].items():
            IS2_atl03_mds[gtx]['geophys_corr'][key] = val[:]

        #-- Getting attributes of included variables
        if ATTRIBUTES:
            #-- Getting attributes of IS2_atl03_mds beam variables
            IS2_atl03_attrs[gtx] = {}
            IS2_atl03_attrs[gtx]['heights'] = {}
            IS2_atl03_attrs[gtx]['geolocation'] = {}
    #             IS2_atl03_attrs[gtx]['bckgrd_atlas'] = {}
            IS2_atl03_attrs[gtx]['geophys_corr'] = {}
            IS2_atl03_attrs[gtx]['Atlas_impulse_response'] = {}
            #-- Global Group Attributes
            for att_name,att_val in fileID[gtx].attrs.items():
                IS2_atl03_attrs[gtx][att_name] = att_val
            #-- ICESat-2 Measurement Group
            for key,val in fileID[gtx]['heights'].items():
                IS2_atl03_attrs[gtx]['heights'][key] = {}
                for att_name,att_val in val.attrs.items():
                    IS2_atl03_attrs[gtx]['heights'][key][att_name]=att_val
            #-- ICESat-2 Geolocation Group
            for key,val in fileID[gtx]['geolocation'].items():
                IS2_atl03_attrs[gtx]['geolocation'][key] = {}
                for att_name,att_val in val.attrs.items():
                    IS2_atl03_attrs[gtx]['geolocation'][key][att_name]=att_val
    #             #-- ICESat-2 Background Photon Rate Group
    #             for key,val in fileID[gtx]['bckgrd_atlas'].items():
    #                 IS2_atl03_attrs[gtx]['bckgrd_atlas'][key] = {}
    #                 for att_name,att_val in val.attrs.items():
    #                     IS2_atl03_attrs[gtx]['bckgrd_atlas'][key][att_name]=att_val
            #-- ICESat-2 Geophysical Corrections Group
            for key,val in fileID[gtx]['geophys_corr'].items():
                IS2_atl03_attrs[gtx]['geophys_corr'][key] = {}
                for att_name,att_val in val.attrs.items():
                    IS2_atl03_attrs[gtx]['geophys_corr'][key][att_name]=att_val

    #-- ICESat-2 spacecraft orientation at time
    IS2_atl03_mds['orbit_info'] = {}
    IS2_atl03_attrs['orbit_info'] = {} if ATTRIBUTES else None
    for key,val in fileID['orbit_info'].items():
        IS2_atl03_mds['orbit_info'][key] = val[:]
        #-- Getting attributes of group and included variables
        if ATTRIBUTES:
            #-- Global Group Attributes
            for att_name,att_val in fileID['orbit_info'].attrs.items():
                IS2_atl03_attrs['orbit_info'][att_name] = att_val
            #-- Variable Attributes
            IS2_atl03_attrs['orbit_info'][key] = {}
            for att_name,att_val in val.attrs.items():
                IS2_atl03_attrs['orbit_info'][key][att_name] = att_val

    #-- number of GPS seconds between the GPS epoch (1980-01-06T00:00:00Z UTC)
    #-- and ATLAS Standard Data Product (SDP) epoch (2018-01-01:T00:00:00Z UTC)
    #-- Add this value to delta time parameters to compute full gps_seconds
    #-- could alternatively use the Julian day of the ATLAS SDP epoch: 2458119.5
    #-- and add leap seconds since 2018-01-01:T00:00:00Z UTC (ATLAS SDP epoch)
    IS2_atl03_mds['ancillary_data'] = {}
    IS2_atl03_attrs['ancillary_data'] = {} if ATTRIBUTES else None
    for key in ['atlas_sdp_gps_epoch']:
        #-- get each HDF5 variable
        IS2_atl03_mds['ancillary_data'][key] = fileID['ancillary_data'][key][:]
        #-- Getting attributes of group and included variables
        if ATTRIBUTES:
            #-- Variable Attributes
            IS2_atl03_attrs['ancillary_data'][key] = {}
            for att_name,att_val in fileID['ancillary_data'][key].attrs.items():
                IS2_atl03_attrs['ancillary_data'][key][att_name] = att_val

    #-- get ATLAS impulse response variables for the transmitter echo path (TEP)
    tep1,tep2 = ('atlas_impulse_response','tep_histogram')
    IS2_atl03_mds[tep1] = {}
    IS2_atl03_attrs[tep1] = {} if ATTRIBUTES else None
    for pce in ['pce1_spot1','pce2_spot3']:
        IS2_atl03_mds[tep1][pce] = {tep2:{}}
        IS2_atl03_attrs[tep1][pce] = {tep2:{}} if ATTRIBUTES else None
        #-- for each TEP variable
        for key,val in fileID[tep1][pce][tep2].items():
            IS2_atl03_mds[tep1][pce][tep2][key] = val[:]
            #-- Getting attributes of included variables
            if ATTRIBUTES:
                #-- Global Group Attributes
                for att_name,att_val in fileID[tep1][pce][tep2].attrs.items():
                    IS2_atl03_attrs[tep1][pce][tep2][att_name] = att_val
                #-- Variable Attributes
                IS2_atl03_attrs[tep1][pce][tep2][key] = {}
                for att_name,att_val in val.attrs.items():
                    IS2_atl03_attrs[tep1][pce][tep2][key][att_name] = att_val

    #-- Global File Attributes
    if ATTRIBUTES:
        for att_name,att_val in fileID.attrs.items():
            IS2_atl03_attrs[att_name] = att_val

    #-- Closing the HDF5 file
    fileID.close()
    #-- Return the datasets and variables
    return (IS2_atl03_mds,IS2_atl03_attrs,IS2_atl03_beams)

# -----------------------------------------------------------------------------------------------

def get_ATL03_x_atc(IS2_atl03_mds, IS2_atl03_attrs, IS2_atl03_beams):
    """
    # Adapted from a notebook by Tyler Sutterly 6/14/2910
    """
    # calculate the along-track and across-track coordinates for ATL03 photons

    Segment_ID = {}
    Segment_Index_begin = {}
    Segment_PE_count = {}
    Segment_Distance = {}
    Segment_Length = {}

    #-- background photon rate
    background_rate = {}

    #-- for each input beam within the file
    for gtx in sorted(IS2_atl03_beams):
        #-- data and attributes for beam gtx
        val = IS2_atl03_mds[gtx]
        val['heights']['x_atc']=np.zeros_like(val['heights']['h_ph'])+np.NaN
        val['heights']['y_atc']=np.zeros_like(val['heights']['h_ph'])+np.NaN
        attrs = IS2_atl03_attrs[gtx]
        #-- ATL03 Segment ID
        Segment_ID[gtx] = val['geolocation']['segment_id']
        n_seg = len(Segment_ID[gtx])
        #-- first photon in the segment (convert to 0-based indexing)
        Segment_Index_begin[gtx] = val['geolocation']['ph_index_beg'] - 1
        #-- number of photon events in the segment
        Segment_PE_count[gtx] = val['geolocation']['segment_ph_cnt']
        #-- along-track distance for each ATL03 segment
        Segment_Distance[gtx] = val['geolocation']['segment_dist_x']
        #-- along-track length for each ATL03 segment
        Segment_Length[gtx] = val['geolocation']['segment_length']
        #-- Transmit time of the reference photon
        delta_time = val['geolocation']['delta_time']

        #-- iterate over ATL03 segments to calculate 40m means
        #-- in ATL03 1-based indexing: invalid == 0
        #-- here in 0-based indexing: invalid == -1
        segment_indices, = np.nonzero((Segment_Index_begin[gtx][:-1] >= 0) &
            (Segment_Index_begin[gtx][1:] >= 0))
        for j in segment_indices:
            #-- index for segment j
            idx = Segment_Index_begin[gtx][j]
            #-- number of photons in segment (use 2 ATL03 segments)
            c1 = np.copy(Segment_PE_count[gtx][j])
            c2 = np.copy(Segment_PE_count[gtx][j+1])
            cnt = c1 + c2
            #-- time of each Photon event (PE)
            segment_delta_times = val['heights']['delta_time'][idx:idx+cnt]
            #-- Photon event lat/lon and elevation (WGS84)
            segment_heights = val['heights']['h_ph'][idx:idx+cnt]
            segment_lats = val['heights']['lat_ph'][idx:idx+cnt]
            segment_lons = val['heights']['lon_ph'][idx:idx+cnt]
            #-- Along-track and Across-track distances
            distance_along_X = np.copy(val['heights']['dist_ph_along'][idx:idx+cnt])
            distance_along_X[:c1] += Segment_Distance[gtx][j]
            distance_along_X[c1:] += Segment_Distance[gtx][j+1]
            distance_along_Y = np.copy(val['heights']['dist_ph_across'][idx:idx+cnt])
            val['heights']['x_atc'][idx:idx+cnt]=distance_along_X
            val['heights']['y_atc'][idx:idx+cnt]=distance_along_Y
