#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import os,re,h5py,pyproj
import glob


def load_data_by_rgt(rgt, path_to_data, product, format = 'hdf5',
                    cycles = ['03','04','05','06','07'], beams = ['gt1l','gt1r','gt2l','gt2r','gt3l','gt3r']):
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

    variables = ['x_atc','times','segment_ids','min_seg_ids','h_li','latitude','longitude','x','y']

    ### extract data from all available cycles
    D_out = {}
    for var in variables:
        D_out[var] = {}

    cycles_this_rgt = []
    for cycle in cycles: # loop over all available cycles
        Di = {}
        for var in variables:
            D_out[var][cycle] = {}

        if not format == 'hdf5':
            raise ValueError('Can only load hdf5 for now.')
        else:
            filenames = glob.glob(os.path.join(path_to_data, f'*{product}_*_{rgt}{cycle}*_003*.h5'))
            error_count=0
            for filename in filenames: # try and load any available files; hopefully is just one
                try:
                    for beam in beams:
                        Di[filename]=atl06_to_dict(filename,'/'+ beam, index=None, epsg=3031, format = 'hdf5')

                        # Times
                        D_out['times'][cycle][beam] = Di[filename]['data_start_utc']
                        # segment ids:
                        D_out['segment_ids'][cycle][beam] = Di[filename]['segment_id']
                        D_out['min_seg_ids'][cycle][beam] = Di[filename]['segment_id'][0]

                        # extract h_li and x_atc, and lat/lons for that section
                        for var in variables:
                            if var in ['times','min_seg_ids','segment_ids']:
                                continue
                            D_out[var][cycle][beam] = Di[filename][var]

                    cycles_this_rgt+=[cycle]
                except KeyError as e:
                    print(f'file {filename} encountered error {e}')
                    error_count += 1

    print('Loaded data from cycles: ' + ','.join(cycles_this_rgt),'for rgt:',rgt)

    return D_out


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

