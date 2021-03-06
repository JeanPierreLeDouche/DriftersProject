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
import scipy.stats as sc

r_e = 6371 * 1e3  # m

data = pickle.load(open("C://Users/Ruben/Documents/CLPH/buoydata.p",'rb'))

# Test to plot a single buoy path

west_border = 310
east_border = 350
north_border = 60
south_border = 45


def lat_corr_distance(lon_distance, lat):
    corrected = np.cos(lat * 180 / np.pi) * lon_distance
    return corrected  # in m


def angle_between_vecs(vector1, vector2):
    x1 = vector1[0]
    y1 = vector1[1]

    x2 = vector2[0]
    y2 = vector2[1]

    dot = x1 * x2 + y1 * y2  # dot product between [x1, y1] and [x2, y2]
    det = x1 * y2 - y1 * x2  # determinant
    angle = np.arctan2(det, dot) * 180 / np.pi  # angle in degrees [-180, 180]

    return angle


buoy_IDs = np.unique(data["ID"])
# dates = np.unique(data["Date"])  # about 16 000

# initialize list where all the angles will go and one temp array
M_angles = pd.DataFrame()
angles_col_np = []

t1 = perf_counter()

# loop through each buoy and then go by time (second for loop)
for ID in buoy_IDs:
    current_buoy_data = data.loc[data['ID'] == ID]
    buoy_lats_lons = current_buoy_data[['Lat', 'Lon']]

    x_previous = np.array([0, 0])
    x_current = np.array([0, 0])
    x_next = np.array([0, 0])

    for time in np.arange(buoy_lats_lons.shape[0] - 2):
        # use three spatial positions at a time to calculate the vectors from
        # first to second and from second to third, giving two consecutive
        # vectors

        x_prev = buoy_lats_lons.iloc[time + 1, :]
        x_curr = buoy_lats_lons.iloc[time, :]
        x_next = buoy_lats_lons.iloc[time + 2, :]

        dx_1 = lat_corr_distance((x_prev[1] - x_curr[1]), x_curr[0])
        dx_2 = lat_corr_distance((x_next[1] - x_curr[1]), x_curr[0])

        dy_1 = (x_prev[0] - x_curr[0])
        dy_2 = (x_next[0] - x_curr[0])

        vec1 = np.array([dx_1, dy_1])
        vec2 = np.asarray([dx_2, dy_2])

        angle = angle_between_vecs(vec1, vec2)

        # angles are stored in a numpy array
        angles_col_np.append(angle)

    # after looping through all timesteps for one buoy the numpy array is copied
    # as a pandas series which can then be appended into a pandas DataFrame.
    # the resulting dataframe has buoy IDs as headers and angles per belonging
    # to this buoy per timestep in each column

    # angles_col_pd = pd.Series(angles_col_np)
    # angles_col_pd = angles_col_pd.mask(angles_col_pd == 0)
    # angles_col_pd = angles_col_pd.dropna(axis=0)
    #
    # M_angles[str(ID)] = angles_col_pd




t2 = perf_counter()
print("calculation took: ", (t2 - t1), " seconds")

# after masking all the trailing zeros we get lists of timestep angles per buoy
# M_angles = M_angles.mask(M_angles == 0)
# M_angles.to_csv()

###--------------------------------------------------------------------------
# %%

# buoy_1 = data.loc[data['ID'] == 34471]
#
# ## plot of 1 buoy that Ruben made
# fig = plt.figure()
#
# ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
# ax.scatter(buoy_1['Lon'], buoy_1['Lat'], s=2)
# ax.coastlines()
# ax.gridlines(draw_labels=True, dms=True)
# ax.set_extent([west_border, east_border, north_border, south_border])
# ax.add_feature(cfeature.OCEAN)
# ax.add_feature(cfeature.LAND)
# ax.add_feature(cfeature.BORDERS)
# ax.add_feature(cfeature.COASTLINE)
# plt.show()
#
# # buoy_1_angles = M_angles
# # # M_his = M_angles
# # # M_his = M_his.dropna(axis = 0)
# # # fig2 = plt.figure()
# # # plt.hist(M_his, 360)
# # # plt.show()
#


x = np.arange(-180,181,1)
a_fit, b_fit = sc.laplace.fit(angles_col_np)
print(a_fit,b_fit)
y = sc.laplace(scale = b_fit)


weights = np.ones_like(angles_col_np) / len(angles_col_np)


fig,(ax1, ax2) = plt.subplots(1,2, sharey=True)
ax2.plot(x, y.pdf(x))
ax2.set_xlabel('Angle (Degrees)')
ax2.set_ylabel('Normalized frequency')


ax1.hist(angles_col_np, 360, weights=weights)
ax1.set_xlabel('Angle (Degrees)')
ax1.set_ylabel('Normalized frequency')

plt.show()

#
# # TO DO: zijn hoeken normaal verdeeld ?
# # paper lezen:
# # decorrelatie tijdsschaal
# # hoe kun je deze hoeken gebruiken om iets over diffusie te zeggen ?




