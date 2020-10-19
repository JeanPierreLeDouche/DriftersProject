import pandas as pd
import numpy as np
import pickle
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from time import perf_counter
import numpy.ma as ma

r_e = 6371 * 1e3 #m 

data = pickle.load(open("BuoyDatabase.p", "rb"))

# print(data.head())

#Test to plot a single buoy path

west_border = 310
east_border = 350
north_border = 60
south_border = 45

def lat_corr_distance(lon_distance, lat):
    corrected =  np.cos(lat * 180 / np.pi )* lon_distance
    return corrected # in m 

def angle_between_vecs(vector1, vector2):
    x1 = vector1[0]
    y1 = vector1[1]
    
    x2 = vector2[0]
    y2 = vector2[1]    
        
    dot = x1*x2 + y1*y2      # dot product between [x1, y1] and [x2, y2]
    det = x1*y2 - y1*x2      # determinant
    angle = np.arctan2(det, dot) * 180 / np.pi  # angle in degrees [-180, 180]
    
    return angle

buoy_IDs = np.unique(data["ID"]) # as of current there are about 10 000
dates = np.unique(data["Date"]) # about 16 000

# initialize list where all the angles will go and one temp array
# M_angles = pd.DataFrame()
# M_angles = np.ndarray(buoy_IDs.shape[0], dtype=object)
M_angles = np.zeros((10000, 16000 ))

t1 = perf_counter()
ID_counter = 0

# loop through each buoy and then go by time (second for loop)
for ID in buoy_IDs[:100]:
    current_buoy_data = data.loc[data['ID'] == ID]
    buoy_lats_lons = current_buoy_data[['Lat','Lon']]
    
    x_previous = np.array([0,0])
    x_current = np.array([0,0])
    x_next = np.array([0,0])
    
    angles_col = np.zeros((16000))

    for time in np.arange(buoy_lats_lons.shape[0]-2):
        
        # use three spatial positions at a time to calculate the vectors from
        # first to second and from second to third, giving two consecutive 
        # vectors 
        
        x_prev = buoy_lats_lons.iloc[time+1,:]
        x_curr = buoy_lats_lons.iloc[time,:]
        x_next = buoy_lats_lons.iloc[time+2,:]
        
        dx_1 = lat_corr_distance((x_prev[1] - x_curr[1]), x_curr[0])
        dx_2 = lat_corr_distance((x_next[1] - x_curr[1]), x_curr[0])
        
        dy_1 =( x_prev[0] - x_curr[0])
        dy_2 =( x_next[0] - x_curr[0]) 
        
        vec1 = np.array([dx_1, dy_1])
        vec2 = np.asarray([dx_2, dy_2])
        
        angle = angle_between_vecs(vec1, vec2)
        
        # angles are stored in a numpy array         
        angles_col[time] = angle
        
    angles_col = np.trim_zeros(angles_col)
    
    # after looping through all timesteps for one buoy the numpy array is copied
    # as a pandas series which can then be appended into a pandas DataFrame. 
    # the resulting dataframe has buoy IDs as headers and angles per belonging
    # to this buoy per timestep in each column
    
    M_angles[ID_counter,:angles_col.shape[0]] = angles_col
    ID_counter = ID_counter + 1

t2 = perf_counter()
print("calculation took: ", (t2-t1), " seconds")

t3 = perf_counter()

zeros_count = 0
for i in range(M_angles.shape[0]):
    for j in range(M_angles.shape[1]):
        if M_angles[i,j] == 0:
            M_angles[i,j] = None
            zeros_count = zeros_count + 1
            
print("counted: ", zeros_count, "zeros")

t4= perf_counter() 
print("dealing with the zeros took:", (t4 - t3), " seconds")

# after masking all the trailing zeros we get lists of timestep angles per buoy    
# M_angles = M_angles.mask(M_angles == 0)    
# M_angles.to_csv()

###--------------------------------------------------------------------------
#%%

# buoy_1 = data.loc[data['ID'] == 34471]

# ## plot of 1 buoy that Ruben made
# fig = plt.figure()

# ax = fig.add_subplot(1,1,1, projection=ccrs.PlateCarree())
# ax.scatter(buoy_1['Lon'],buoy_1['Lat'], s=2)
# ax.coastlines()
# ax.gridlines(draw_labels=True, dms=True)
# ax.set_extent([west_border, east_border, north_border, south_border])
# ax.add_feature(cfeature.OCEAN)
# ax.add_feature(cfeature.LAND)
# ax.add_feature(cfeature.BORDERS)
# ax.add_feature(cfeature.COASTLINE)
# plt.show()

#%%

plt.hist(M_angles, 360)


=======
#Data before filtering
print(data.shape)

#Filtering data

data = data.loc[data['Lat'] > south_border]
data = data.loc[data['Lat'] < north_border]
data = data.loc[data['Lon'] > west_border]
data = data.loc[data['Lon'] < east_border]

#Data after filtering
print(data.shape)

#Checking the amount of different buoys in the area
print(f'The amount of buoys in the area is: {len(set(data["ID"]))}')

# buoy_1 = data.loc[data['ID'] == 72615]
#
# fig = plt.figure()
#
# ax = fig.add_subplot(1,1,1, projection=ccrs.PlateCarree())
# ax.scatter(buoy_1['Lon'],buoy_1['Lat'], s=2)
# ax.coastlines()
# ax.gridlines(draw_labels=True, dms=True)
# ax.set_extent([west_border, east_border, north_border, south_border])
# ax.add_feature(cfeature.OCEAN)
# ax.add_feature(cfeature.LAND)
# ax.add_feature(cfeature.BORDERS)
# ax.add_feature(cfeature.COASTLINE)
#
#
# plt.show()
