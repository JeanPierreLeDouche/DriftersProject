
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

# path = r'C:/Users/Ruben/Documents/CLPH/MAIO/DrifterProject/buoydata_15001_mar20.dat'
# data = pd.read_table(path, header=None, sep='\s+')

west_border = 310
east_border = 350
north_border = 60
south_border = 45

print('Data loaded, started filtering...')

data.columns = ["ID", "Month", "Day", "Year", "Lat", "Lon", "SST(Deg C)", "VE(CM/S)", "VN(CM/S)",
                "SPD(CM/S)", "VAR.LAT", "VAR.LON", "VAR.TEMP"]

to_drop = ["SST(Deg C)","VE(CM/S)","VN(CM/S)", "VAR.LAT", "VAR.LON", "VAR.TEMP"]

data.drop(to_drop, inplace = True, axis = 'columns')

data.loc[data['Lat'] > 90, 'Lat'] = np.nan
data.loc[data['Lon'] > 360, 'Lon'] = np.nan

data = data.loc[data['Year'] > 2011]

### introduce single date and hourly time format

data["Hour"] = (data["Day"]-np.floor(data["Day"]) )*24 # convert decimals from "Day" colunn into hours
data["Date"] = pd.to_datetime(data[["Year", "Month", "Day"]], yearfirst =True) + pd.to_timedelta(data["Hour"], unit = 'h')

# remove the now unneccesary columns
data.drop(["Year", "Month", "Day", "Hour"], inplace=True, axis = 'columns')

# sort by date to get a timeseries
data = data.sort_values(by="Date")

# This file is saved in the same folder as the Python file, an absolute path can be added to save at another place

data_filtered = data.loc[data['Lat'] > south_border]
data_filtered = data_filtered.loc[data['Lat'] < north_border]
data_filtered = data_filtered.loc[data['Lon'] > west_border]
data_filtered = data_filtered.loc[data['Lon'] < east_border]

buoy_IDs = np.unique(data["ID"])
buoy_IDs_f = np.unique(data_filtered["ID"])

good_ids = []

for i in enumerate(buoy_IDs_f):
    if len(data.loc[data['ID'] == i[1]]) == len(data_filtered.loc[data['ID'] == i[1]]):
        good_ids.append(i[1])

data_correct = data_filtered[data_filtered.ID.isin(good_ids)]


print('Data filtered and sorted, now saving as pickle file...')

# This file is saved in the same folder as the Python file, an absolute path can be added to save at another place

pickle.dump(data_correct, open(r"C:\Users\Gebruiker\Documents\Climate_Physics\Year2\MAIO\Driftersproject\buoydata.p","wb"))


