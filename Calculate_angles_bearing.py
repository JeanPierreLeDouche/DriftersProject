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
from scipy.stats import vonmises


# n = 100  # number of buoys to be calculated

data = pickle.load(
    open(r'C:\Users\Gebruiker\Documents\Climate_Physics\Year2\MAIO\Driftersproject\ALL_buoydata.p', "rb"))

def bearing_angle(p1, p2):
    Y = np.sin((p2[0] - p1[0])*(np.pi/180)) * np.cos(p2[1]*(np.pi/180))
    X = np.cos(p1[1]*(np.pi/180))*np.sin(p2[1]*(np.pi/180)) - np.sin(p1[1]*(np.pi/180)) * np.cos(p2[1]*(np.pi/180)) * np.cos((p2[0]-p1[0])*(np.pi/180))

    b = np.arctan2(X,Y)

    return b

def tau(delta, psi):
    # from Visser 2008
    if psi.any() == 0:
        value = ''
    else:
        value = delta / (1 - psi)
    return value


# as defined in Visser 2008
def diffusion(v, delta, psi, n=2):
    # from Visser 2008
    if psi.any() == 0:
        value = ''
    else:
        value = (1. / n) * ((v ** 2 * delta) / (1 - psi))
    return value

def normalize(list1):
    arr = np.asarray(list1)
    mean = np.mean(arr)
    std = np.std(arr)
    arr = (arr - mean)/std    
    return arr


buoy_IDs = np.unique(data["ID"])

t1 = perf_counter()

# create two grids that are 180 lists that contain 360 lists that contain a
# single list. Using np.array is not suited because this asks for fixed
# "depth" of the grid which is variable in this application

L_pos = []
L_ang = []
L_ang_NH = []
L_ang_SH = []

lat_points = 180
lon_points = 361

angles_grid = [[[] for x in range(lon_points)] for x in range(lat_points + 1)]
speeds_grid = [[[] for x in range(lon_points)] for x in range(lat_points + 1)]

# loop through each buoy and then go by time (second for loop)
for ID in enumerate(buoy_IDs):
    print(ID[0])
    current_buoy_data = data.loc[data['ID'] == ID[1]]
    buoy_lons_lats = current_buoy_data[['Lon', 'Lat']]
    buoy_speed = current_buoy_data[['SPD(CM/S)']]
    x_previous = np.array([0, 0])
    x_current = np.array([0, 0])
    x_next = np.array([0, 0])
    bear = []
    angles = []
    angles_col = np.zeros((current_buoy_data.shape[0]))

    for time in np.arange(buoy_lons_lats.shape[0] - 2):
        # first and last speed are 999.999, we exclude the first by indexing
        # the speed with time+1 which runs until the total time -2 thereby
        # also excluding the last value

        v_curr = buoy_speed.iloc[time + 1]

        # use three spatial positions at a time to calculate the vectors from
        # first to second and from second to third, giving two consecutive
        # vectors

        x_prev = buoy_lons_lats.iloc[time, :] #0th element lon, 1th element lan
        x_curr = buoy_lons_lats.iloc[time + 1, :]
        x_next = buoy_lons_lats.iloc[time + 2, :]


        point_1 = [x_prev[0],x_prev[1]]
        point_2 = [x_curr[0], x_curr[1]]
        point_3 = [x_next[0], x_next[1]]

        bearing1 = bearing_angle(point_1, point_2)
        bearing2 = bearing_angle(point_2, point_3)
        # bearings = np.array([bearing1, bearing2])

        bear.append(bearing1)
        
        angle = (bearing2 - bearing1) * (180 / np.pi)  # in degrees
        
        if angle > 180:
            # print('True')
            angle -= 360
        elif angle < - 180:
            angle += 360

        angles.append(angle)
        
        if x_curr[1] > 0:
            L_ang_NH.append(angle)
        elif x_curr[1] < 0: 
            L_ang_SH.append(angle)    
        
        # lats and lons are rounded down to their nearest integer
        lon = int(x_prev[0])
        lat = int(x_prev[1])

        new_lat = lat + 90
        new_lon = lon

        # autocor_grid[new_lat][new_lon].append(autoc)
        angles_grid[new_lat][new_lon].append(angle)
        speeds_grid[new_lat][new_lon].append(v_curr)

        L_ang.append(angle)

    #Comparing possibilities for psi
    # print(bear)
    # print(f'The autocorrelation is: {sm.tsa.acf(angles, nlags=1)[1]}')
    # angles = np.asarray(angles)
    # angles = np.cos(angles*(np.pi/180))
    # angles = np.mean(angles)

t2 = perf_counter()
print(" angles calculation took: ", (t2 - t1), " seconds")

# %%

psi_grid = np.zeros((180, 360))
avg_speed_grid = np.zeros((180, 360))

for lat in range(0, 180):
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
        avg_speed_grid[lat][lon] = avg_speed / 100

tau_grid = tau(6 * 3600, psi_grid)
D_grid = diffusion(avg_speed_grid, 6 * 3600, psi_grid)




### vonmises fitting and plots
#%%
def misesfitplot(data, title):
    # calculate statistical metrics
    mean =np.round( np.mean(data), 2)
    std = np.round(np.std(data), 2)
    datapoints = len(data)
    
    #normalize data
    norm_data = normalize(data) - 0.03
    
    # fit von mises distr.
    misfit =  vonmises.fit(norm_data, fscale=1)
    kappa = misfit[0]
    loc = misfit[1]
    scale = misfit[2]

    x = np.linspace(vonmises.ppf(0.001, kappa),
                    vonmises.ppf(0.999, kappa), 100)
    
    # actual plotting
    plt.title(title, fontsize =20)
    plt.hist(norm_data, 360, density=True, label='normalized data')
    plt.text(-3, 0.7, f'mean: {mean} $\degree$', fontsize = 18)
    plt.text(-3, 0.65, f'std: {std} $\degree$', fontsize = 18)
    plt.text(-3, 0.60, f'datapoints: {datapoints}', fontsize = 18)
    plt.plot(x, vonmises.pdf(x, kappa=kappa, loc=loc, scale=scale), color='r', lw = 3)
    plt.ylabel('(normalized) frequency', fontsize = 15)
    plt.xlabel('(normalized) angle', fontsize = 15)
    plt.yticks(fontsize =15)
    plt.xticks(fontsize =15)
    
    plt.show()
    return 
#%%
misesfitplot(np.asarray(L_ang)*-1, "Histogram of all angles")
misesfitplot(np.asarray(L_ang_SH) *-1, 'Histogram of SH angles')
misesfitplot(np.asarray(L_ang_NH) * -1, 'Histogram of NH angles')



# x = []
# y = []

###Single buoy plot and angles
#
# for i in range(len(L_pos)):
#     x.append(L_pos[i][0])
#     y.append(L_pos[i][1])
#
# x = x[:9]
# y = y[:9]
#
# # ax1 = plt.axes(projection=ccrs.PlateCarree())
# # ax1.coastlines()
# # ax1.plot(x, y)
# # ax1.text(x[0], y[0], 'start')
# # ax1.scatter(x[0], y[0], color='green')
# # ax1.text(x[len(x) - 1], y[len(y) - 1], 'end')
# # ax1.scatter(x[len(x) - 1], y[len(y) - 1], color='r')
# # ax1.set_title('Single buoy path')
# # for i in range(1, len(x) - 1):
# #     ax1.text(x[i], y[i], str(i))
# # # ax1.set_ylim(36.8, 37.6)
# # # ax1.set_xlim(161.5, 162)
# #
# # plt.show()
# #
# # ax2 = plt.scatter(range(1,len(L_ang[:9])-1), (np.asarray(L_ang[:7])), ls='-.', lw=3)
# # plt.title('Counter-clockwise angles')
# # plt.ylabel(r'Angle ($\degree$)')
# # plt.xlabel('step #')
# #
# # plt.show()
# # %%

###_____________________________________________________________________________________________________________________
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
# %%
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
#
#%%
# fit normal distributions
# x = np.arange(-180,181,1)
# norm_a, norm_b = sc.norm.fit(np.asarray(L_ang)) # see if we can make this better
# y_norm = sc.norm(norm_a, norm_b)

# plt.figure()
# plt.ylabel('frequency')
# plt.xlabel(r'Angle ($\theta$)')
# plt.hist(np.asarray(L_ang), 360, density = True, stacked =True) # check if normalization is done correctly
# plt.plot(x, y_norm.pdf(x), color = 'r')
# # plt.show()

# plt.plot(x, y.pdf(x), color = 'r')

# plt.legend(['Fitted normal distribution', 'Data'])

# plt.show()

# s, p = sc.normaltest(L_ang)
# print(p)


# #%%_________________________________________________________________________________________________

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