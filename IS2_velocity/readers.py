#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import os,re,h5py,pyproj
import glob
import pandas as pd

def atl06_to_dict(filename, beam, field_dict=None, index=None, epsg=None, format = 'hdf5'):
    """
    Read selected datasets from an ATL06 file

    Input arguments:
        filename: ATl06 file to read
        beam: a string specifying which beam is to be read (ex: gt1l, gt1r, gt2l, etc)
        field_dict: A dictinary describing the fields to be read
                keys give the group names to be read,
                entries are lists of datasets within the groups
        index: which entries in each field to read
        epsg: an EPSG code specifying a projection (see www.epsg.org).  Good choices are:
            for Greenland, 3413 (polar stereographic projection, with Greenland along the Y axis)
            for Antarctica, 3031 (polar stereographic projection, centered on the Pouth Pole)
        format: what format the files are in. Default = hd5f. No other formats currently supported.
    Output argument:
        D6: dictionary containing ATL06 data.  Each dataset in
            dataset_dict has its own entry in D6.  Each dataset
            in D6 contains a numpy array containing the
            data
    """
    if format == 'hdf5':
        if field_dict is None:
            field_dict={None:['latitude','longitude','h_li', 'atl06_quality_summary'],\
                        'ground_track':['x_atc','y_atc'],\
                        'fit_statistics':['dh_fit_dx', 'dh_fit_dy']}
        D={}
        file_re=re.compile('ATL06_(?P<date>\d+)_(?P<rgt>\d\d\d\d)(?P<cycle>\d\d)(?P<region>\d\d)_(?P<release>\d\d\d)_(?P<version>\d\d).h5')
        with h5py.File(filename,'r') as h5f:
            for key in field_dict:
                for ds in field_dict[key]:
                    if key is not None:
                        ds_name=beam+'/land_ice_segments/'+key+'/'+ds
                    else:
                        ds_name=beam+'/land_ice_segments/'+ds
                    if index is not None:
                        D[ds]=np.array(h5f[ds_name][index])
                    else:
                        D[ds]=np.array(h5f[ds_name])
                    if '_FillValue' in h5f[ds_name].attrs:
                        bad_vals=D[ds]==h5f[ds_name].attrs['_FillValue']
                        D[ds]=D[ds].astype(float)
                        D[ds][bad_vals]=np.NaN
            D['data_start_utc'] = h5f['/ancillary_data/data_start_utc'][:]
            D['delta_time'] = h5f['/' + beam + '/land_ice_segments/delta_time'][:]
            D['segment_id'] = h5f['/' + beam + '/land_ice_segments/segment_id'][:]
        if epsg is not None:
            xy=np.array(pyproj.proj.Proj(epsg)(D['longitude'], D['latitude']))
            D['x']=xy[0,:].reshape(D['latitude'].shape)
            D['y']=xy[1,:].reshape(D['latitude'].shape)
        temp=file_re.search(filename)
        D['rgt']=int(temp['rgt'])
        D['cycle']=int(temp['cycle'])
        D['beam']=beam
    return D



### Some functions
# MISSING HERE: mask by data quality?
def load_data_by_rgt(rgt, path_to_data, product, filter_type = 'running_average', running_avg_window = None, format = 'hdf5'):
    """
    rgt: repeat ground track number of desired data
    filter_type: Name of desired filter type.
        'running_average' = a centered running avergae filter of smoothing_window_size will be used
        None = return raw data, no filter
    running_average_window: Default None
    dx: desired spacing. NOT USED, removed from function (should be retrieved from data itself)
    path_to_data:
    product: ex., ATL06
    format, ex 'hdf5'
    """

    # hard code these for now:
    cycles = ['03','04','05','06','07'] # not doing 1 and 2, because don't overlap exactly
    beams = ['gt1l','gt1r','gt2l','gt2r','gt3l','gt3r']

    ### extract data from all available cycles
    x_atc = {}
    lats = {}
    lons = {}
    h_li_raw = {} # unsmoothed data; equally spaced x_atc, still has nans
    h_li_raw_NoNans = {} # unsmoothed data; equally spaced x_atc, nans filled with noise
    h_li = {} # smoothed data, equally spaced x_atc, nans filled with noise
    h_li_diff = {}
    times = {}
    min_seg_ids = {}
    segment_ids = {}
    x_ps= {}
    y_ps= {}

    cycles_this_rgt = []
    for cycle in cycles: # loop over all available cycles
        Di = {}
        x_atc[cycle] = {}
        lats[cycle] = {}
        lons[cycle] = {}
        h_li_raw[cycle] = {}
        h_li_raw_NoNans[cycle] = {}
        h_li[cycle] = {}
        h_li_diff[cycle] = {}
        times[cycle] = {}
        min_seg_ids[cycle] = {}
        segment_ids[cycle] = {}
        x_ps[cycle]= {}
        y_ps[cycle]= {}

        if format == 'hdf5':

            filenames = glob.glob(os.path.join(path_to_data, f'*{product}_*_{rgt}{cycle}*_003*.h5'))
            error_count=0


            for filename in filenames: # try and load any available files; hopefully is just one
                try:
                    for beam in beams:
                        Di[filename]=atl06_to_dict(filename,'/'+ beam, index=None, epsg=3031, format = 'hdf5')



                        times[cycle][beam] = Di[filename]['data_start_utc']

                        # extract h_li and x_atc, and lat/lons for that section
                        x_atc_tmp = Di[filename]['x_atc']
                        h_li_tmp = Di[filename]['h_li']#[ixs]
                        lats_tmp = Di[filename]['latitude']
                        lons_tmp = Di[filename]['longitude']
                        x_ps_tmp = Di[filename]['x']
                        y_ps_tmp= Di[filename]['y']


                        # segment ids:
                        seg_ids = Di[filename]['segment_id']
                        min_seg_ids[cycle][beam] = seg_ids[0]
                        #print(len(seg_ids), len(x_atc_tmp))

                        # make a monotonically increasing x vector
                        # assumes dx = 20 exactly, so be carefull referencing back
                        ind = seg_ids - np.nanmin(seg_ids) # indices starting at zero, using the segment_id field, so any skipped segment will be kept in correct location
                        # TASK: Change dx to depend on data loaded, not be hard coded as 20 (as in line below)
                        x_full = np.arange(np.max(ind)+1) * 20 + x_atc_tmp[0]
                        h_full = np.zeros(np.max(ind)+1) + np.NaN
                        h_full[ind] = h_li_tmp
                        lats_full = np.zeros(np.shape(x_full)) * np.nan
                        lats_full[ind] = lats_tmp
                        lons_full = np.zeros(np.shape(x_full)) * np.nan
                        lons_full[ind] = lons_tmp
                        x_ps_full = np.zeros(np.shape(x_full)) * np.nan
                        x_ps_full[ind] = x_ps_tmp
                        y_ps_full = np.zeros(np.shape(x_full)) * np.nan
                        y_ps_full[ind] = y_ps_tmp

                        ## save the segment id's themselves, with gaps filled in
                        segment_ids[cycle][beam] = np.zeros(np.max(ind)+1) + np.NaN
                        segment_ids[cycle][beam][ind] = seg_ids


                        x_atc[cycle][beam] = x_full
                        h_li_raw[cycle][beam] = h_full # preserves nan values
                        lons[cycle][beam] = lons_full
                        lats[cycle][beam] = lats_full
                        x_ps[cycle][beam] = x_ps_full
                        y_ps[cycle][beam] = y_ps_full

                        ### fill in nans with noise h_li datasets
                #                         h = ma.array(h_full,mask =np.isnan(h_full)) # created a masked array, mask is where the nans are
                #                         h_full_filled = h.mask * (np.random.randn(*h.shape)) # fill in all the nans with random noise

                        ### interpolate nans in pandas
                        # put in dataframe for just this step; eventually rewrite to use only dataframes?
                        data = {'x_full': x_full, 'h_full': h_full}
                        df = pd.DataFrame(data, columns = ['x_full','h_full'])
                        #df.plot(x='x_full',y='h_full')
                        # linear interpolation for now
                        df['h_full'].interpolate(method = 'linear', inplace = True)
                        h_full_interp = df['h_full'].values
                        h_li_raw_NoNans[cycle][beam] = h_full_interp # has filled nan values

                        ### Apply any filters; differentiate
                        dx = np.median(x_full[1:] - x_full[:-1])
                        if filter_type == 'running_average':
                            if running_avg_window == None:
                                running_avg_window = 100

                            h_filt = filt(x1 = x_full, h1 = h_full_interp, dx = dx, filter_type = 'running_average', running_avg_window = running_avg_window)
                            h_li[cycle][beam] = h_filt

                            # differentiate that section of data
                            h_diff = differentiate(x_full, h_filt) #(h_filt[1:] - h_filt[0:-1]) / (x_full[1:] - x_full[0:-1])
                        elif filter_type == None:
                            h_li[cycle][beam] = h_full_interp
                            h_diff = differentiate(x_full, h_full_interp) # (h_full_interp[1:] - h_full_interp[0:-1]) / (x_full[1:] - x_full[0:-1])

                        h_li_diff[cycle][beam] = h_diff



                        #print(len(x_full), len(h_full), len(lats_full), len(seg_ids), len(h_full_interp), len(h_diff))


                    cycles_this_rgt+=[cycle]
                except KeyError as e:
                    print(f'file {filename} encountered error {e}')
                    error_count += 1

    print('Cycles available: ' + ','.join(cycles_this_rgt))
    return x_atc, lats, lons, h_li_raw, h_li_raw_NoNans, h_li, h_li_diff, \
            times, min_seg_ids, segment_ids, cycles_this_rgt, x_ps, y_ps




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
