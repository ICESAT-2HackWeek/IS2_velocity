# from icepyx import icesat2data as ipd
import os, glob, re, h5py, sys, pyproj
import matplotlib as plt
import shutil
import numpy as np
from pprint import pprint
from astropy.time import Time
from scipy.signal import correlate, detrend
import pandas as pd
import matplotlib.pyplot as plt
# %matplotlib widget
import pointCollection as pc




##### Where is the data: ################
### Where are the data to be processed
#datapath = '/home/jovyan/shared/surface_velocity/FIS_ATL06'

#local data path
# datapath = '/media/rag110/ADATA SD700/ICESat2/download/FIS'

map_data_root = '/Users/grace/Dropbox/Cornell/projects/003/FIS_data/'
datapath = '/Users/grace/Dropbox/Cornell/projects/003/ATL06_0848/data/'
# datapath = #'/home/jovyan/shared/surface_velocity/FIS_ATL06'
ATL06_files=glob.glob(os.path.join(datapath, '*.h5'))


### Where to save the results
# out_path = '/home/jovyan/shared/surface_velocity/ATL06_out2/'
out_path = '/Users/grace/Dropbox/Cornell/projects/003/ATL06_0848/out/'

#local out_path Rodrigo
# out_path = '/media/rag110/ADATA SD700/ICESat2/output/FIS/'







