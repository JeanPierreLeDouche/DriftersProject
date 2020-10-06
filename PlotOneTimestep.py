# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 13:21:31 2020

@author: Gebruiker
"""

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from sys import argv
import pandas as pd
import numpy as np
import pickle 

#filepath = r"C:\Users\Gebruiker\Documents\Climate_Physics\Year2\MAIO\Driftersproject\buoydata_15001_mar20"
#filepath1 = 'buoydata_15001_mar20.dat'
#
#data = pd.read_table(filepath1,header = None, sep= '\s+') 
#
#headernames = ["ID","Month" , "Day", "Year", "Lat","Lon","SST(Deg C)","VE(CM/S)","VN(CM/S)",
#                 "SPD(CM/S)","VAR.LAT","VAR.LON","VAR.TEMP"]

#data.columns = headernames

data = pickle.load(open("BuoyDatabase.p", "rb"))

### fin unique ID's for buoy's
buoy_IDs = np.unique(data["ID"])

# get location data for one timestep

loc_data = data[["ID","Month", "Day", "Year", "Lat", "Lon" ]]

loc_data = loc_data.sort_values(by=["Year", "Month", "Day"]) # sort data on date 

hourdata = (loc_data["Day"]-np.floor(loc_data["Day"]) )*24 # convert decimals from days into hours
loc_data["Hour"] = hourdata 

dates = pd.to_datetime(loc_data[["Year", "Month", "Day"]], yearfirst =True) + pd.to_timedelta(loc_data["Hour"], unit = 'h')
    
loc_data.drop(["Year", "Month", "Day", "Hour"], inplace=True, axis = 'columns')
loc_data["Date"] = dates



#### plotting all locations 
#
#ax = plt.axes(projection=ccrs.PlateCarree())
#ax.coastlines()
#
#plt.show()

    