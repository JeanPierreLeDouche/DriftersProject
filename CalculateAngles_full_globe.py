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
data = pickle.load(open(r"C:\Users\Gebruiker\Documents\Climate_Physics\Year2\MAIO\Driftersproject\ALL_buoydata.p", "rb"))

def lat_corr_distance(lon_distance, lat):    
    # circles of equal latitude become smaller towards the poles therefore we need
    # to correct for this when calculating absolute distances from distances between
    # longitudes
    corrected =  np.cos(lat * 180 / np.pi )* lon_distance
    return corrected # in m 

def angle_between_vecs(vector1, vector2):
    x1 = vector1[0]
    y1 = vector1[1]
    
    x2 = vector2[0]
    y2 = vector2[1]    
        
    dot = x1*x2 + y1*y2      # dot product between [x1, y1] and [x2, y2]
    det = x1*y2 - y1*x2      # determinant
    angle = np.arctan2(det, dot) #angle in radians
    
    return angle

def tau(delta, psi):
    # from Visser 2008
    if psi.any() == 0:
        value = ''
    else:
        value = delta/(1-psi)   
    return value

# as defined in Visser 2008
def diffusion(v, delta, psi, n=2 ):
   # from Visser 2008
    if psi.any() == 0:
        value = ''
    else: 
        value = (1./n) * ((v**2 * delta) /(1 - psi) )
    return value


buoy_IDs = np.unique(data["ID"]) 

t1 = perf_counter()

L_angles = []
L_speeds = []

# create two grids that are 180 lists that contain 360 lists that contain a 
# single list. Using np.array is not suited because this asks for fixed
# "depth" of the grid which is variable in this application

lat_points = 180
lon_points = 360
  
angles_grid = [[[] for x in range(lon_points)] for x in range(lat_points+1)] 
speeds_grid = [[[] for x in range(lon_points)] for x in range(lat_points+1)] 

# loop through each buoy and then go by time (second for loop)
for ID in enumerate(buoy_IDs[:1000]):
    current_buoy_data = data.loc[data['ID'] == ID[1]]  
    buoy_lons_lats = current_buoy_data[['Lon', 'Lat']]
    buoy_speed = current_buoy_data[['SPD(CM/S)']]
    print(ID[1])
    x_previous = np.array([0,0])
    x_current = np.array([0,0])
    x_next = np.array([0,0])
    
    angles_col = np.zeros((current_buoy_data.shape[0]))

    for time in np.arange(buoy_lons_lats.shape[0]-2):
        
        # first and last speed are 999.999, we exclude the first by indexing 
        # the speed with time+1 which runs until the total time -2 thereby
        # also excluding the last value 
        
        v_curr = buoy_speed.iloc[time+1]
        
        # use three spatial positions at a time to calculate the vectors from
        # first to second and from second to third, giving two consecutive 
        # vectors 
        
        x_prev = buoy_lons_lats.iloc[time+1,:]
        x_curr = buoy_lons_lats.iloc[time,:]
        x_next = buoy_lons_lats.iloc[time+2,:]
        
        dx_1 = lat_corr_distance((x_prev[0] - x_curr[0]), x_curr[1])
        dx_2 = lat_corr_distance((x_next[0] - x_curr[0]), x_curr[1])
        
        dy_1 = (x_prev[1] - x_curr[1])
        dy_2 = (x_next[1] - x_curr[1]) 
        
        vec1 = np.array([dx_1, dy_1])
        vec2 = np.array([dx_2, dy_2])
        
        angle = angle_between_vecs(vec1, vec2) # in radians
        
        # lats and lons are rounded down to their nearest integer
        lon = int(x_prev[0]) 
        lat = int(x_prev[1])
        
        # new_lat = (lat - 90)*-1
        # new_lon = (lon+180)%360
        
        new_lat = lat + 90
        new_lon = lon
        
        angles_grid[new_lat][new_lon].append(angle)
        speeds_grid[new_lat][new_lon].append(v_curr)
    
    # diagnostics 
    
    if ID[0]==1000:
        t11 = perf_counter()
        print('1000 buoys calculated after: ', (t11-t1), ' seconds' )
    if ID[0] == 2500:
        t12 = perf_counter()
        print('2500 buoys calculated after: ', (t12-t1), ' seconds')
    if ID[0] == 5000:
        t13 = perf_counter()
        print('5000 buoys calculated after: ', (t13-t1), ' seconds')
    if ID[0] == 7500: 
        t14 = perf_counter()
        print('7500 buoys calculated after: ', (t14-t1), 'seconds')
        
t2 = perf_counter()
print(" angles calculation took: ", (t2-t1), " seconds")

#%%

# calculating psi and average speeds at every gridpoint

psi_grid = np.zeros((180,360))
avg_speed_grid = np.zeros((180,360))

for lat in range(0 , 180):
    for lon in range(0, 360): 

        angle_grid_point_slice = np.asarray(angles_grid[lat][lon])
        if angle_grid_point_slice.any() == True:
            phi = np.mean(np.cos(angle_grid_point_slice))
        else: 
            phi = np.NaN
        psi_grid[lat][lon] = phi

        speed_grid_point_slice = np.asarray(speeds_grid[lat][lon])
        if speed_grid_point_slice.any() == True:
            avg_speed = np.mean(speed_grid_point_slice)
        else: 
            avg_speed = np.NaN
        avg_speed_grid[lat][lon] = avg_speed/100

# using the psi and avg. v values for the full grid we calculate tau and D
# for the full grid
        
tau_grid = tau(6*3600, psi_grid)
D_grid = diffusion(avg_speed_grid, 6 * 3600, psi_grid)

#%%
# Box with sea of South Africa, for reproducing Rühs 2018
west_border = 10
east_border = 50
north_border = -19
south_border = -49

lon = np.linspace(0,360 , 360)
lat = np.linspace(-90, 90, 180)

data_crs = ccrs.PlateCarree()

ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_title(r'$\psi$ ')
ax.coastlines()
# ax.legend()
ax.set_global()
ax.gridlines(crs=ccrs.PlateCarree(), linewidth=0.5, color='black', 
    draw_labels=True, alpha=0.5, linestyle='--')
ax.contourf(lon, lat, psi_grid, transform=data_crs, color = 'red')


#%%

plt.figure()
plt.title('Correlation timescale')
plt.contourf(tau_grid)
plt.colorbar()
plt.show()

# plt.figure()
# plt.title('Diffusivity')
# plt.contourf(D_grid)
# plt.colorbar()
# plt.show()

###---------------------------------------------------------------------------
" below here should be all post processing stuff that is omitted for speed"
# # as defined in Visser 2008
# def tau(delta, psi):
#     return delta/(1-psi)

# # as defined in Visser 2008
# def diffusion(v, delta, psi, n=2 ):
#     return (1./n) * ((v**2 * delta) /(1 - psi) )
#%%
# # calculate total psi (cos of angle between vectors from before) 
# ensemble_psi = np.mean(np.asarray(np.cos(L_angles)))

# # calculate psi per buoy
# buoys_psi_list = []
# for i in enumerate(buoy_IDs):
#     buoys_psi_list.append(np.mean(np.cos(np.trim_zeros(np.asarray(M_angles[i[0],:])))))

# # convert speeds from cm/s to m/s    
# speeds = np.mean(np.asarray(L_speeds) /100)
# speeds = speeds.flatten()

# # correlation timescale tau, for whole dataset and individual buoys
# ensemble_tau = tau(6 * 3600, ensemble_psi)
# buoys_tau = tau(6 * 3600, np.asarray(buoys_psi_list))

# # diffusion D, for whole dataset and individual buoys
# ensemble_D = diffusion(np.mean(speeds), 6 * 3600, ensemble_psi)
# buoys_D = diffusion(speeds, 6 * 3600, np.asarray(buoys_psi_list))   

# print('Average tau and D are: ', ensemble_tau,' ',ensemble_D, 'respectively')
# print('avg speed is : ', np.mean(speeds), " m/s")
# print('For longitudes between', west_border, 'and', east_border)



# ### ---------------------------------------------------------------------------

# # fit normal and laplace distributions
# x = np.arange(-180,181,1)
# a_fit, b_fit = sc.laplace.fit(np.asarray(L_angles) * 180/np.pi)
# norm_a, norm_b = sc.norm.fit(np.asarray(L_angles) *180/np.pi) # see if we can make this better
# y = sc.laplace(scale = b_fit)
# y_norm = sc.norm(norm_a, norm_b)
# # weights = np.ones_like(np.asarray(L_angles) / len(L_angles))

# #%%

# plt.figure()
# plt.ylabel('frequency')
# plt.xlabel(r'Angle ($\theta$)')
# plt.hist(np.asarray(L_angles) * 180/np.pi, 360, density = True, stacked =True) # check if normalization is done correctly


# plt.plot(x, y.pdf(x), color = 'r')
# plt.plot(x, y_norm.pdf(x), color = 'y')
# plt.legend(['Fitted Laplace distribution','Fitted normal distribution', 'Data'])

# plt.show()


# #%%

# plt.figure()
# plt.title('Correlation timescale')
# plt.hist(np.log10(buoys_tau), 50)
# plt.hist(np.log10(ensemble_tau), 50, color = 'r')
# plt.ylabel('frequency')
# plt.xlabel(r'log ($\tau$)')  #units!
# plt.show()

# #%%

# plt.figure()
# plt.title('Diffusivity')
# plt.hist(np.log10(buoys_D), bins=50)
# plt.hist(np.log10(ensemble_D),15, color = 'r')
# plt.ylabel('frequency')
# plt.xlabel(r'log (D)')
# plt.show()

