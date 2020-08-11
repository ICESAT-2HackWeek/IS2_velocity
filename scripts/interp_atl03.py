#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
# Import some of the scripts that we have written
import glob, sys
sys.path.append("/Users/benhills/Software/MyGitRepos/surface_velocity/scripts")
sys.path.append("/Users/benhills/Software/MyGitRepos/surface_velocity/readers")
from read_HDF5_ATL03 import read_HDF5_ATL03
from get_ATL03_x_atc import get_ATL03_x_atc

import warnings
warnings.filterwarnings('ignore')

atl3_dir = '/Users/benhills/Google Drive File Stream/Shared drives/IceSat2_Surface_Velocity/Shared_Data/FIS_0848_ATL03/ATL03_files/'
files = glob.glob(atl3_dir+'**/ATL03*.h5')


for file_name in files[2:]:

    # read the IS2 data with Tyler's ATL03 reader:
    IS2_atl03_mds, IS2_atl03_attrs, IS2_atl03_beams = read_HDF5_ATL03(file_name)
    # add x_atc to the ATL03 data structure (this function adds to the LS2_ATL03_mds dictionary)
    get_ATL03_x_atc(IS2_atl03_mds, IS2_atl03_attrs, IS2_atl03_beams)

    beams = ['gt1l', 'gt1r', 'gt2l', 'gt2r', 'gt3l', 'gt3r']

    beam = beams[0]

    #-- select the 2l beam from ATL03
    D3 = IS2_atl03_mds[beam]
    LMH=D3['heights']['signal_conf_ph'][:,3] >= 2
    atl03_along = D3['heights']['x_atc'][LMH]
    atl03_height = D3['heights']['h_ph'][LMH]

    xstep = 5.
    xwin = 10.
    xs = np.arange(np.nanmin(atl03_along),np.nanmax(atl03_along),xstep)
    hs = np.empty_like(xs)
    print('Total Length:',len(xs))

    for i,x in enumerate(xs):
        if i%100 == 0:
            print(i)
        idxs = np.logical_and(atl03_along>x-xwin,atl03_along<x+xwin)
        xfit = atl03_along[idxs]
        hfit = atl03_height[idxs]
        idx = ~np.isnan(xfit) & ~np.isnan(hfit)
        if ~np.any(idx):
            hs[i] = np.nan
        else:
            fit = np.polyfit(xfit,hfit,1);
            hs[i] = np.polyval(fit,x);

    out_name = file_name[:-3] + '_' + beams[0] + '_spacing_' + str(round(xstep)) + '_window_' + str(round(xwin))
    np.save(out_name,np.array([xs,hs]))
