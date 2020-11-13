# Import the basic libraries
import numpy as np
import matplotlib.pyplot as plt
# %matplotlib notebook
import os,re,h5py,pyproj, glob

# Import the basic libraries
import numpy as np
import matplotlib.pyplot as plt

# This cell: put the versions of codes I'm working on first in path
import sys
print(sys.path)
working_module_directory = '/Users/grace/Dropbox/Cornell/projects/003/git_repo'
if not sys.path[0] == working_module_directory:
    sys.path.insert(0, working_module_directory)
print(sys.path)

# As an example, import a function from the ICESat-2 surface velocity library
from IS2_velocity.correlation_processing import calculate_velocities

from IS2_velocity.readers_atl06 import load_data_by_rgt

########## testinga bunch of data
# where is data
data_path = '/Users/grace/Dropbox/Cornell/projects/003/FIS_03_04/raw_ATL06/'
ATL06_files=glob.glob(os.path.join(data_path, '*.h5'))


# where to save raw data
save_path = '/Users/grace/Dropbox/Cornell/projects/003/FIS_03_04/corr_results/'

# List of available rgts
rgts = {}
rgts_list = []
for filepath in ATL06_files:
    filename = filepath.split('/')[-1]
    rgt = filename.split('_')[3][0:4]
    track = filename.split('_')[3][4:6]
#     print(rgt,track)
    if not rgt in rgts.keys():
        rgts[rgt] = []
        rgts[rgt].append(track)
    else:
        rgts[rgt].append(track)
    rgts_list.append(rgt)


# all rgt values in our study are are in rgts.keys()
print(rgts.keys())

# available tracks for each rgt are in rgts[rgt]; ex.:
print(rgts['0848'])

# Calculate velocity between cycles 3 and 4
cycle1 = '03'
cycle2 = '04'
beams = ['gt1l', 'gt1r', 'gt2l', 'gt2r', 'gt3l', 'gt3r']

# Loop over rgts; Load data; preprocess; run cross correlations; write out
from IS2_velocity.preprocessing import *
from IS2_velocity.saving import *
from IS2_velocity.saving import *
from IS2_velocity.plotting import plot_measures_along_track_comparison


# saving options
write_out_path = save_path
write_out_prefix = 'run20201112'
map_data_root = '/Users/grace/Dropbox/Cornell/projects/003/FIS_data/' # where map data are stored; different on each computer
spatial_extent = np.array([-65, -86, -55, -81])
product = 'ATL06'
dx = 20
segment_length = 2000
search_width = 1000
along_track_step = 100
max_percent_nans = 10

error_rgts = []
for ir, rgt in enumerate(rgts_list):
    print(rgt)
    if cycle1 in rgts[rgt] and cycle2 in rgts[rgt]:
        # if there is data for both cycles (there is not necessarily)

        try:
                     # load data
            data = load_data_by_rgt(rgt, data_path, product='ATL06', format='hdf5')

            # preprocess
            data = fill_seg_ids(data)
            data = interpolate_nans(data)
            data = filt(data, filter_type='running_average', running_avg_window=60)
            data = differentiate(data)

            # run correlations
            data = calculate_velocities(data, cycle1, cycle2, beams)

            # save results
            save_is2_velocity(data, rgt, write_out_path, write_out_prefix, product, map_data_root,
                              cycle1, cycle2, dx, segment_length, search_width, along_track_step,
                              max_percent_nans, spatial_extent)
            #
            # plotting = True
            # if plotting:
            #     plot_out_location = write_out_path
            #     correlation_threshold = 0.65
            #     velocity_number = 0
            #     spatial_extent = np.array([-65, -86, -55, -81])
            #     epsg = 3031
            #     plot_measures_along_track_comparison(rgt, beams, write_out_path, correlation_threshold, spatial_extent,
            #                                          plot_out_location, map_data_root, velocity_number, epsg)
        except Exception:
            print('problem with data from rgt ' + str(rgt))
            error_rgts += [rgt]
            print(Exception)


processed_files = glob.glob(write_out_path + '/*hdf5')
for ir, rgt in enumerate(rgts_list):
    print(rgt)
    if cycle1 in rgts[rgt] and cycle2 in rgts[rgt]:
        try:
            plot_out_location = write_out_path
            correlation_threshold = 0.65
            velocity_number = 0
            spatial_extent = np.array([-65, -86, -55, -81])
            epsg = 3031
            plot_measures_along_track_comparison(rgt, beams, write_out_path, correlation_threshold, spatial_extent,
                                                 plot_out_location, map_data_root, velocity_number, epsg)
        except:
            pass

# problem rgts: 0513, 0842, 0491, 0772, 0771, 0878, 0720, 0482, 0833, 0492

##################################





# path to data, relative to folder /notebooks
data_dir = 'data/'#'../data/'
rgt = '0848'

# Load data; This step loads raw data, interpolates to constant spacing, filters if requested, and
# differentiates
filter_type = 'running_average'
running_avg_window = 100

data = load_data_by_rgt(rgt, data_dir, product = 'ATL06',format = 'hdf5')

# preprocessing
from IS2_velocity.preprocessing import *

data = fill_seg_ids(data)
data = interpolate_nans(data)
data = filt(data)
data = differentiate(data)

## calc veloces
from IS2_velocity.correlation_processing import calculate_velocities

# Calculate velocity between cycles 3 and 4
cycle1 = '03'
cycle2 = '04'
beams = ['gt1l','gt1r','gt2l','gt2r','gt3l','gt3r']

data = calculate_velocities(data, cycle1, cycle2, beams)

################ OLD
# As an example, import a function from the ICESat-2 surface velocity library
# from IS2_velocity.correlation_processing import velocity
# help(velocity)

##### TEST READER atl06
# Import the reader script
# from IS2_velocity.readers import atl06_to_dict

# read in dictionaries from two different cycles
data_dir = 'data/' # '../data/'
fn_1 = 'processed_ATL06_20190822153035_08480411_003_01.h5'
D1=atl06_to_dict(data_dir+fn_1,'/gt2l', index=None, epsg=3031, format='hdf5')
fn_2 = 'processed_ATL06_20190523195046_08480311_003_01.h5'
D2=atl06_to_dict(data_dir+fn_2,'/gt2l', index=None, epsg=3031, format='hdf5')

# Plot the landice elevation along the pass.
plt.figure(figsize=(8,4))

plt.plot(D1['x_atc']/1000.,D1['h_li'],c='indianred')
plt.plot(D2['x_atc']/1000.,D2['h_li'],c='steelblue')
plt.ylabel('Elevation (m)')
plt.xlabel('Along-Track Distance (km)')

plt.tight_layout()

# Get segment ids from the loaded dictionaries
x1,h1 = fill_seg_ids(D1['x_atc'],D1['h_li'],D1['segment_id'])
x2,h2 = fill_seg_ids(D2['x_atc'],D2['h_li'],D2['segment_id'])

# Smooth and differentiate the elevation product (this is a preprocessing step)
running_avg_window = 100
h1_filt = filt(x1, h1, running_avg_window, filter_type = 'running_average')
h2_filt = filt(x2, h2, running_avg_window, filter_type = 'running_average')
dh1 = differentiate(x1,h1)
dh2 = differentiate(x2,h2)
dh1_filt = differentiate(x1,h1_filt)
dh2_filt = differentiate(x2,h2_filt)



plt.figure(figsize=(8,6))

# Plot smoothed surface elevation
plt.subplot(211)
plt.tick_params(labelbottom=False,bottom=False)
plt.plot(x1/1000.,h1,c='grey')
plt.plot(x1/1000.,h1_filt,c='k')
plt.ylabel('Elevation (m)')

# Plot the surface Slope
plt.subplot(212)
plt.plot(x1/1000.,dh1,c='grey')
plt.plot(x1/1000.,dh1_filt,c='k')
plt.xlabel('Along-Track Distance (km)')
plt.ylabel('Surface Slope (m/m)')

plt.tight_layout()

# Calculate time offset
dt = time_diff(D1,D2)

### Control the correlation step:
segment_length = 2000 # meters, how wide is the window we are correlating in each step
search_width = 1000 # meters, how far in front of and behind the window to check for correlation
along_track_step = 100 # meters; how much to jump between each consecutivevelocity determination
max_percent_nans = 10 # Maximum % of segment length that can be nans and still do the correlation step

x_atc, lats, lons, h_li_raw, h_li_raw_NoNans, h_li, h_li_diff, times, min_seg_ids, segment_ids, cycles_this_rgt, x_ps, y_ps = \
    load_data_by_rgt(rgt = '0848', path_to_data = data_dir, product = 'ATL06', \
                     filter_type = 'running_average', running_avg_window = running_avg_window, \
                     format = 'hdf5')


cycles_this_rgt = sorted(cycles_this_rgt)
cycle1 = cycles_this_rgt[0]
cycle2 = cycles_this_rgt[1]

# dt = time_diff # Currently coded to use hdf5 data dictionary;
### Timing of each cycle in the current velocity determination
t1_string = times[cycle1]['gt1l'][0].astype(str)  # figure out later if just picking hte first one it ok
t1 = Time(t1_string)
t2_string = times[cycle2]['gt1l'][0].astype(str)  # figure out later if just picking hte first one it ok
t2 = Time(t2_string)

### Which beams to process
beams = ['gt1l','gt1r','gt2l','gt2r','gt3l','gt3r']

### Elapsed time between cycles
dt = (t2 - t1).jd  # difference in julian days

###
rgt = '0848'

dx = 20
velocities, correlations, lags, midpoints_x_atc, midpoints_xy, midpoints_lons, midpoints_lats = \
    velocity(x_atc, h_li_diff, lats, lons, segment_ids, beams, cycle1, cycle2, \
             segment_length, search_width, along_track_step, max_percent_nans, dx)



from matplotlib.gridspec import GridSpec

beam = 'gt1l'
x1 = x_atc['03'][beam]
x2 = x_atc['04'][beam]
h1 = h_li['03'][beam]
h2 = h_li['04'][beam]
dh1 = h_li_diff['03'][beam]
dh2 = h_li_diff['04'][beam]
vel_xs = midpoints_x_atc[rgt][beam]
velocs = velocities[rgt][beam]

plt.figure(figsize=(8,4))
gs = GridSpec(2,2)

# Plot the elevation profiles again
plt.subplot(gs[0,0])
plt.tick_params(bottom=False,labelbottom=False)
plt.plot(x1/1000.-29000,h1,'.',c='indianred')
plt.plot(x2/1000.-29000,h2,'.',c='steelblue',ms=3)
plt.ylabel('Elevation (m)')
plt.title('ATL06',fontweight='bold')
plt.xlim(80,580)

# Plot the slopes again
plt.subplot(gs[1,0])
plt.tick_params(bottom=False,labelbottom=False)
plt.plot(x1/1000.-29000,dh1,'.',c='indianred')
plt.plot(x2/1000.-29000,dh2,'.',c='steelblue',ms=3)
plt.ylim(-.05,.05)
plt.ylabel('Surface Slope (m/m)')
plt.xlim(80,580)

# Plot the calculated velocities along track
ax5 = plt.subplot(gs[0,1])
plt.plot(vel_xs/1000.-29000,velocs,'.',c='k',label='ATL06')
plt.ylabel('Velocity (m/yr)')
plt.xlabel('Along-Track Distance (km)')
plt.xlim(80,580)
plt.ylim(-500,1500)

plt.tight_layout()








