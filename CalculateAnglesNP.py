# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 09:42:15 2020

@author: Gebruiker
"""

import pandas as pd
import numpy as np
import pickle
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from time import perf_counter

data = pickle.load(open("BuoyDatabase.p", "rb"))

#Test to plot a single buoy path

west_border = 310
east_border = 350
north_border = 60
south_border = 45

buoy_IDs = np.unique(data["ID"]) # as of current there are about 10 000
dates = np.unique(data["Date"]) # about 16 000

# initialize list where all the angles will go and one temp array
M_angles = pd.DataFrame()
angles_col_np = np.zeros((16000))

# loop through each buoy and then go by time (second for loop)
t_start = perf_counter()

for ID in buoy_IDs[:0]:
    current_buoy_data = data.loc[data['ID'] == ID]
    buoy_lats_lons = current_buoy_data[['Lat','Lon']]
     
    x_previous = pd.DataFrame([[0,0]])
    x_current = pd.DataFrame([[0],[0]])
            
    
    for time in np.arange(buoy_lats_lons.shape[0]-2):
        
        # use three spatial positions at a time to calculate the vectors from
        # first to second and from second to third, giving two consecutive 
        # vectors 
        t_0 = perf_counter()
        x_previous = buoy_lats_lons.iloc[time,:]
        x_current = buoy_lats_lons.iloc[time+1,:]
        x_next = buoy_lats_lons.iloc[time+2,:]
        
        vec1 = np.asarray(x_current - x_previous)
        vec2 = np.asarray(x_next - x_current)
            
        # using that the dot product of a,b is defined : a.b = |a||b|cos(angle)
        # where the angle is the angle between the vectors in this case, we can
        # exploit this definition by using that angle = arccos(a.b / |a||b|)
        
        dot = np.dot(vec2, vec1)
        
        normalized_dot = dot/ (np.linalg.norm(vec1)*np.linalg.norm(vec2))
        
        angle = np.arccos(normalized_dot) * 180 /np.pi 
        #  
        
        # angles are stored in a numpy array         
        angles_col_np[time] = angle
        

    # after looping through all timesteps for one buoy the numpy array is copied
    # as a pandas series which can then be appended into a pandas DataFrame. 
    # the resulting dataframe has buoy IDs as headers and angles per belonging
    # to this buoy per timestep in each column
    
    # NOTE: there should be a more elegent way of doing this 
    
    angles_col_pd = pd.Series(angles_col_np)
    M_angles[str(ID)] = angles_col_pd 



# after masking all the trailing zeros we get lists of timestep angles per buoy    
M_angles = M_angles.mask(M_angles == 0)    
M_angles.to_csv()



###--------------------------------------------------------------------------
#buoy_1 = data.loc[data['ID'] == 72615]

### plot of 1 buoy that Ruben made
#fig = plt.figure()
#
#ax = fig.add_subplot(1,1,1, projection=ccrs.PlateCarree())
#ax.scatter(buoy_1['Lon'],buoy_1['Lat'], s=2)
#ax.coastlines()
#ax.gridlines(draw_labels=True, dms=True)
#ax.set_extent([west_border, east_border, north_border, south_border])
#ax.add_feature(cfeature.OCEAN)
#ax.add_feature(cfeature.LAND)
#ax.add_feature(cfeature.BORDERS)
#ax.add_feature(cfeature.COASTLINE)
#
#
#plt.show()
