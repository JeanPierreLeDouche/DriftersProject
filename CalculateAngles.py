import pandas as pd
import numpy as np
import pickle
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from time import perf_counter
import numpy.ma as ma
import scipy.stats as sc

r_e = 6371 * 1e3 #m 

data = pickle.load(open(r"C:\Users\Gebruiker\Documents\Climate_Physics\Year2\MAIO\Driftersproject\buoydata.p", "rb"))

# west_border = 310
# east_border = 350
west_border = 330
east_border = 350
north_border = 60
south_border = 45

data = data.loc[data['Lon'] > west_border]
data = data.loc[data['Lon'] < east_border]

# circles of equal latitude become smaller towards the poles therefore we need
# to correct for this when calculating absolute distances from distances between
# longitudes
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

# as defined in Visser 2008
def tau(delta, psi):
    return delta/(1-psi)

# as defined in Visser 2008
def diffusion(v, delta, psi, n=2 ):
    return (1./n) * ((v**2 * delta) /(1 - psi) )

buoy_IDs = np.unique(data["ID"]) 
dates = np.unique(data["Date"]) 

# initialize array where all the angles will go,  
M_angles = np.zeros((100, 4000 ))

t1 = perf_counter()
ID_counter = 0 

L_angles = []
L_speeds = []

# loop through each buoy and then go by time (second for loop)
for ID in buoy_IDs:
    current_buoy_data = data.loc[data['ID'] == ID]  
    buoy_lats_lons = current_buoy_data[['Lat','Lon']]
    buoy_speed = current_buoy_data[['SPD(CM/S)']]
    
    # in taking the avg speed we exclude the first and last value because these 
    # are 999.999. Assuming this is some kind of measurement issue. 
    # avg_speed = buoy_speed[1:-1].mean()
    
    x_previous = np.array([0,0])
    x_current = np.array([0,0])
    x_next = np.array([0,0])
    
    angles_col = np.zeros((current_buoy_data.shape[0]))

    for time in np.arange(buoy_lats_lons.shape[0]-2):
        
        # use three spatial positions at a time to calculate the vectors from
        # first to second and from second to third, giving two consecutive 
        # vectors 
        
        x_prev = buoy_lats_lons.iloc[time+1,:]
        x_curr = buoy_lats_lons.iloc[time,:]
        x_next = buoy_lats_lons.iloc[time+2,:]
        
        dx_1 = lat_corr_distance((x_prev[1] - x_curr[1]), x_curr[0])
        dx_2 = lat_corr_distance((x_next[1] - x_curr[1]), x_curr[0])
        
        dy_1 = (x_prev[0] - x_curr[0])
        dy_2 = (x_next[0] - x_curr[0]) 
        
        vec1 = np.array([dx_1, dy_1])
        vec2 = np.asarray([dx_2, dy_2])
        
        angle = angle_between_vecs(vec1, vec2) * np.pi/180 # in radians
        
        # angles are stored in a numpy array         
        angles_col[time] = angle
        # also stored in one big list
        L_angles.append(angle)
        
        # save avg speed per buoy  

        L_speeds.append(buoy_speed.iloc[time+1] )
    
    # arrays of angles are put into a big matrix, each row containing angles
    # belonging to a specific buoy
    M_angles[ID_counter,:angles_col.shape[0]] = angles_col

    ID_counter = ID_counter + 1

t2 = perf_counter()
print("calculation took: ", (t2-t1), " seconds")

#%%
# calculate total psi (cos of angle between vectors from before) 
ensemble_psi = np.mean(np.asarray(np.cos(L_angles)))

# calculate psi per buoy
buoys_psi_list = []
for i in enumerate(buoy_IDs):
    buoys_psi_list.append(np.mean(np.cos(np.trim_zeros(np.asarray(M_angles[i[0],:])))))

# convert speeds from cm/s to m/s    
speeds = np.mean(np.asarray(L_speeds) /100)
speeds = speeds.flatten()

# correlation timescale tau, for whole dataset and individual buoys
ensemble_tau = tau(6 * 3600, ensemble_psi)
buoys_tau = tau(6 * 3600, np.asarray(buoys_psi_list))

# diffusion D, for whole dataset and individual buoys
ensemble_D = diffusion(np.mean(speeds), 6 * 3600, ensemble_psi)
buoys_D = diffusion(speeds, 6 * 3600, np.asarray(buoys_psi_list))   

print('Average tau and D are: ', ensemble_tau,' ',ensemble_D, 'respectively')
print('avg speed is : ', np.mean(speeds), " m/s")
print('For longitudes between', west_border, 'and', east_border)

### ---------------------------------------------------------------------------

# fit normal and laplace distributions
x = np.arange(-180,181,1)
a_fit, b_fit = sc.laplace.fit(np.asarray(L_angles) * 180/np.pi)
norm_a, norm_b = sc.norm.fit(np.asarray(L_angles) *180/np.pi)
y = sc.laplace(scale = b_fit)
y_norm = sc.norm(norm_a, norm_b)
# weights = np.ones_like(np.asarray(L_angles) / len(L_angles))

#%%

plt.figure()
plt.ylabel('frequency')
plt.xlabel(r'Angle ($\theta$)')
plt.hist(np.asarray(L_angles) * 180/np.pi, 360, density = True, stacked =True)

plt.plot(x, y.pdf(x), color = 'r')
plt.plot(x, y_norm.pdf(x), color = 'y')
plt.legend(['Fitted Laplace distribution','Fitted normal distribution', 'Data'])

plt.show()


#%%

plt.figure()
plt.title('Correlation timescale')
plt.hist(np.log10(buoys_tau), 50)
plt.hist(np.log10(ensemble_tau), 50, color = 'r')
plt.ylabel('frequency')
plt.xlabel(r'log ($\tau$)')
plt.show()

#%%

plt.figure()
plt.title('Diffusivity')
plt.hist(np.log10(buoys_D), bins=50)
plt.hist(np.log10(ensemble_D),15, color = 'r')
plt.ylabel('frequency')
plt.xlabel(r'log (D)')
plt.show()


