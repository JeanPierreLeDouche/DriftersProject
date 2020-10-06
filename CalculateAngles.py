# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 09:42:15 2020

<<<<<<< Updated upstream
@author: Gebruiker
"""

=======
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
>>>>>>> Stashed changes
