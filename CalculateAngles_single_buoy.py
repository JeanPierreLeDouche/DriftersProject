import pandas as pd
import numpy as np
import pickle
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from time import perf_counter
import numpy.ma as ma
import scipy.stats as sc

n = 100 # number of buoys to be calculated

r_e = 6371 * 1e3 #m 
circ_e = 2 * np.pi * r_e
circ_e_per_deg =  circ_e/360

data = pickle.load(open(r"C:\Users\Gebruiker\Documents\Climate_Physics\Year2\MAIO\Driftersproject\ALL_buoydata.p", "rb"))

def lat_corr_distance(lon_distance, lat):    
    # circles of equal latitude become smaller towards the poles therefore we need
    # to correct for this when calculating absolute distances from distances between
    # longitudes
    corrected =  (r_e * np.cos(lat * np.pi / 180 ) * 2 * np.pi) / 360 * lon_distance
    return corrected # in m
 
def bearing(point1, point2):
    #phi,labda
    
    
    x = np.sin(point2[0] - point1[0]) * np.cos(point2[1])  
    y = np.cos(point1[1]) * np.sin(point2[1]) - np.sin(point1[1]) * np.cos(point2[1]) * np.cos((point2[0]-point1[0]) )
    theta = np.arctan2(y,x)
    bear = (theta * 180 /np.pi + 360) % 360 
    
    return bear
   
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

# create two grids that are 180 lists that contain 360 lists that contain a 
# single list. Using np.array is not suited because this asks for fixed
# "depth" of the grid which is variable in this application

L_pos =[]
L_ang = []

lat_points = 180
lon_points = 360
  
angles_grid = [[[] for x in range(lon_points)] for x in range(lat_points+1)] 
speeds_grid = [[[] for x in range(lon_points)] for x in range(lat_points+1)] 

# loop through each buoy and then go by time (second for loop)
for ID in enumerate(buoy_IDs[1:2]):
    current_buoy_data = data.loc[data['ID'] == ID[1]]  
    buoy_lons_lats = np.asarray(current_buoy_data[['Lon', 'Lat' ]])
    buoy_speed = current_buoy_data[['SPD(CM/S)']]
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
        
        x_prev = buoy_lons_lats[time,:]
        x_curr = buoy_lons_lats[time+1,:]
        x_next = buoy_lons_lats[time+2,:]
        
        bear_1 =  bearing(x_prev, x_curr)
        bear_2 = bearing(x_curr, x_next)

        angle = bear_2 - bear_1 # in radians
        
        # lats and lons are rounded down to their nearest integer
        lon = int(x_prev[0]) 
        lat = int(x_prev[1])
        
        new_lat = lat + 90
        new_lon = lon
        
        angles_grid[new_lat][new_lon].append(angle)
        speeds_grid[new_lat][new_lon].append(v_curr)
        
        L_pos.append(x_prev)
        L_ang.append(angle)
        
t2 = perf_counter()
print(" angles calculation took: ", (t2-t1), " seconds")

#%%

x = []
y = []

for i in range(len(L_pos)):
    x.append(L_pos[i][0])
    y.append(L_pos[i][1])
    
x = x[:11]
y = y[:11]    


ax1 = plt.axes(projection=ccrs.PlateCarree())
ax1.coastlines()
ax1.plot(x, y)
ax1.text(x[0], y[0], 'start')
ax1.scatter(x[0], y[0], color = 'green')
ax1.text(x[len(x)-1], y[len(y)-1], 'end')
ax1.scatter(x[len(x)-1], y[len(y)-1], color = 'r')

for i in range(1,len(x)-2):
    ax1.text(x[i+1], y[i+1], str( np.around(L_ang[i] , 0)))

ax1.set_ylim(36.8,37.6)
ax1.set_xlim(161.5, 162)
plt.show()



#%%

# plt.figure()
# plt.title('Correlation timescale')
# plt.contourf(tau_grid)
# plt.colorbar()
# plt.show()

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

