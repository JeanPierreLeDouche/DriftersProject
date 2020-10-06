# -*- coding: utf-8 -*-


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


buoy_1 = data.loc[data['ID'] == 72615]


fig = plt.figure()

ax = fig.add_subplot(1,1,1, projection=ccrs.PlateCarree())
ax.scatter(buoy_1['Lon'],buoy_1['Lat'], s=2)
ax.coastlines()
ax.gridlines(draw_labels=True, dms=True)
ax.set_extent([west_border, east_border, north_border, south_border])
ax.add_feature(cfeature.OCEAN)
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.BORDERS)
ax.add_feature(cfeature.COASTLINE)


plt.show()


