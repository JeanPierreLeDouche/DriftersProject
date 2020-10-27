import pandas as pd
import numpy as np
import pickle
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from time import perf_counter
import numpy.ma as ma
import scipy.stats as sc
import statsmodels.api as sm
from scipy.stats import norm
import scipy

data = pickle.load(
    open(r'C:\Users\Gebruiker\Documents\Climate_Physics\Year2\MAIO\Driftersproject\ALL_buoydataforRuhs.p', "rb"))


D_grid = [[[] for x in range(lon_points + 1)] for x in range(lat_points + 1)]


### Plotting

# fig, ax  = plt.subplots(1,1,figsize=(10,5), subplot_kw={'projection': ccrs.PlateCarree()})
#
# cbar_ax = fig.add_axes([0.9, 0.2, 0.03, 0.7])
# fig.subplots_adjust(hspace=0, wspace=0, top=0.90, right=0.8)
#
# bounds_c = [1e2, 5e2, 1e3, 5e3, 1e4, 5e4, 1e5]
#
# ax.set_title(r'D ')
# ax.coastlines()
# # ax.legend()
# ax.set_global()
# ax.gridlines(crs=ccrs.PlateCarree(), linewidth=0.5, color='black',
#     draw_labels=True, alpha=0.5, linestyle='--')
# c = ax.contourf(lon, lat, D_grid, bounds_c, norm=col.LogNorm(), transform=data_crs, extend = 'both')
#
# plt.colorbar(c, cax=cbar_ax, ax=ax)
#
# plt.savefig('C://Users/Ruben/Documents/CLPH/MAIO/dplotfull')
#
# plt.show()