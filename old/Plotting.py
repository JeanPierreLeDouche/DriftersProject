-*- coding: utf-8 -*-


import pandas as pd
import numpy as np
import pickle
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature




data = pickle.load(open("BuoyDatabase.p", "rb"))

# print(data.head())

#Test to plot a single buoy path

west_border = 310
east_border = 350
north_border = 60
south_border = 45

#Data before filtering
print(f'The raw data consists of:{data.shape} data points')

#Filtering data

data = data.loc[data['Lat'] > south_border]
data = data.loc[data['Lat'] < north_border]
data = data.loc[data['Lon'] > west_border]
data = data.loc[data['Lon'] < east_border]


#Data after filtering
print(f'The filtered data consists of:{data.shape} data points')


#Checking the amount of different buoys in the area
buoy_IDs = list(set(data["ID"]))
num_bouys = len(buoy_IDs)
print(f'The amount of buoys in the area is: {num_bouys}')

data_lat = data["Lat"]
data_lon = data["Lon"]


data_crs = ccrs.PlateCarree()

fig = plt.figure()

ax = fig.add_subplot(1,1,1, projection=ccrs.Orthographic())

ax.gridlines(draw_labels=True, dms=True)
ax.set_extent([west_border, east_border, north_border, south_border])
ax.add_feature(cfeature.OCEAN)
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.BORDERS)
ax.add_feature(cfeature.COASTLINE)

ax.scatter(np.array(data["Lon"]), np.array(data["Lat"]), s=0.05, transform = data_crs)

plt.show()