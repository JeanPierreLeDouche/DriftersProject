import pandas as pd
import numpy as np
import pickle


"""
This Python file loads in the buoy database, adds the correct headers, selects only the needed columns and
replaces incorrect data by NaN values. This modified database is then saved as a pickle file, so it can be loaded
quickly for further use.
"""

# Replace this path with the correct path
path = r'C:\Users\Gebruiker\Documents\Climate_Physics\Year2\MAIO\Driftersproject\buoydata_15001_mar20.dat'
data = pd.read_table(path, header=None, sep='\s+')

print('Data loaded, startedfiltering...')

data.columns = ["ID", "Month", "Day", "Year", "Lat", "Lon", "SST(Deg C)", "VE(CM/S)", "VN(CM/S)",
                "SPD(CM/S)", "VAR.LAT", "VAR.LON", "VAR.TEMP"]

to_drop = ["SST(Deg C)", "VAR.LAT", "VAR.LON", "VAR.TEMP"]

data.drop(to_drop, inplace = True, axis = 'columns')

data.loc[data['Lat'] > 90, 'Lat'] = np.nan
data.loc[data['Lon'] > 360, 'Lon'] = np.nan

data = data.loc[data['Year'] > 2011]
data = data.loc[data['Year'] < 2018]

# ### introduce single date and hourly time format

# data["Hour"] = (data["Day"]-np.floor(data["Day"]) )*24 # convert decimals from "Day" colunn into hours
# data["Date"] = pd.to_datetime(data[["Year", "Month", "Day"]], yearfirst =True) + pd.to_timedelta(data["Hour"], unit = 'h')

# # remove the now unneccesary columns
# data.drop(["Year", "Month", "Day", "Hour"], inplace=True, axis = 'columns')

# # sort by date to get a timeseries
# data = data.sort_values(by="Date")

print('Data filtered and sorted, now saving as pickle file...')

pickle.dump(data, open("BuoyDatabaseforRuhs.p","wb"))

