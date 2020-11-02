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
from matplotlib import ticker, cm
import matplotlib.colors as col
from matplotlib import colorbar
# n = 100  # number of buoys to be calculated

data = pickle.load(
    open(r"C:\Users\Ruben\Documents\CLPH\MAIO\ALL_buoydata.p", "rb"))



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


buoy_IDs = np.unique(data["ID"])

t1 = perf_counter()

# create two grids that are 180 lists that contain 360 lists that contain a
# single list. Using np.array is not suited because this asks for fixed
# "depth" of the grid which is variable in this application

L_pos = []
L_ang = []

lat_points = 180
lon_points = 360

angles_grid = [[[] for x in range(lon_points + 1)] for x in range(lat_points + 1)]
speeds_grid = [[[] for x in range(lon_points + 1)] for x in range(lat_points + 1)]

# loop through each buoy and then go by time (second for loop)
loop_size = buoy_IDs[:100]

for ID in enumerate(loop_size):
    if ID[0] % (len(loop_size) // 100) == 0:
        print(f'At {int(100*(ID[0] / len(loop_size)))}%')
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





        # print(angle)


        if angle > 180:
            # print('True')
            angle -= 360
        elif angle < - 180:
            angle += 360

        angle *= -1

        angles.append(angle)

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


    # print(f'The mean of the cosines is: {angles}')
    # diagnostics

    if ID[0] == 1000:
        t11 = perf_counter()
        print('1000 buoys calculated after: ', (t11 - t1), ' seconds')
    if ID[0] == 2500:
        t12 = perf_counter()
        print('2500 buoys calculated after: ', (t12 - t1), ' seconds')
    if ID[0] == 5000:
        t13 = perf_counter()
        print('5000 buoys calculated after: ', (t13 - t1), ' seconds')
    if ID[0] == 7500:
        t14 = perf_counter()
        print('7500 buoys calculated after: ', (t14 - t1), 'seconds')



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

pickle.dump(psi_grid, open("C://users/Ruben/Documents/CLPH/MAIO/fullworldpsi.p","wb"))
pickle.dump(D_grid, open("C://users/Ruben/Documents/CLPH/MAIO/fullworldD.p","wb"))
pickle.dump(tau_grid, open("C://users/Ruben/Documents/CLPH/MAIO/fullworldtau.p","wb"))

x = []
y = []

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
# # plt.title('Clockwise angles')
# # plt.ylabel(r'Angle ($\degree$)')
# # plt.xlabel('step #')
# #
# # plt.show()
# # %%

###_____________________________________________________________________________________________________________________


# ax1 = plt.axes(projection=ccrs.PlateCarree())
# # ax1.coastlines()
# # ax1.plot(x, y)
# # ax1.text(x[0], y[0], 'start')
# # ax1.scatter(x[0], y[0], color='green')
# # ax1.text(x[len(x) - 1], y[len(y) - 1], 'end')
# # ax1.scatter(x[len(x) - 1], y[len(y) - 1], color='r')
# # ax1.set_title('Single buoy path')

# ax1 = plt.axes(projection=ccrs.PlateCarree())
# ax1.coastlines()


# plt.figure()
#
# plt.title('Psi')
# plt.xlabel('Lat')
# plt.ylabel('Lon')
# plt.contourf(psi_grid[:][:])
# plt.colorbar()
#
# plt.show()

# west_border = 10
# east_border = 50
# north_border = -19
# south_border = -49

lon = np.linspace(0,360 , 360)
lat = np.linspace(-90, 90, 180)

data_crs = ccrs.PlateCarree()



fig, ax  = plt.subplots(1,1,figsize=(10,5), subplot_kw={'projection': ccrs.PlateCarree()})

cbar_ax = fig.add_axes([0.9, 0.2, 0.03, 0.7])
fig.subplots_adjust(hspace=0, wspace=0, top=0.90, right=0.8)
a = ax.contourf(lon, lat, psi_grid, transform=data_crs)
ax.set_title(r'$\psi$ ')
ax.coastlines()
# ax.legend()
ax.set_global()

ax.gridlines(crs=ccrs.PlateCarree(), linewidth=0.5, color='black',
    draw_labels=True, alpha=0.5, linestyle='--')


plt.colorbar(a, cax=cbar_ax, ax=ax)

plt.show()

###

fig, ax  = plt.subplots(1,1,figsize=(10,5), subplot_kw={'projection': ccrs.PlateCarree()})

bounds_b = [1e4, 5e4, 1e5, 5e5, 1e6]


cbar_ax = fig.add_axes([0.9, 0.2, 0.03, 0.7])
fig.subplots_adjust(hspace=0, wspace=0, top=0.90, right=0.8)
b = ax.contourf(lon, lat, tau_grid, bounds_b, transform=data_crs, extend = 'both', norm=col.LogNorm())
ax.set_title(r'$\tau$ ')
ax.coastlines()
# ax.legend()
ax.set_global()

ax.gridlines(crs=ccrs.PlateCarree(), linewidth=0.5, color='black',
    draw_labels=True, alpha=0.5, linestyle='--')

plt.colorbar(b, cax=cbar_ax, ax=ax)


plt.show()


# ax = plt.axes(projection=ccrs.PlateCarree())
# ax.set_title(r'$\tau$ ')
# ax.coastlines()
# # ax.legend()
# ax.set_global()
# ax.gridlines(crs=ccrs.PlateCarree(), linewidth=0.5, color='black',
#     draw_labels=True, alpha=0.5, linestyle='--')
# b = ax.contourf(lon, lat, tau_grid, norm=col.LogNorm(), transform=data_crs)
#
# plt.colorbar(b)
#
# plt.show()

fig, ax  = plt.subplots(1,1,figsize=(10,5), subplot_kw={'projection': ccrs.PlateCarree()})

cbar_ax = fig.add_axes([0.9, 0.2, 0.03, 0.7])
fig.subplots_adjust(hspace=0, wspace=0, top=0.90, right=0.8)

bounds_c = [1e2, 5e2, 1e3, 5e3, 1e4, 5e4, 1e5]

ax.set_title(r'D ')
ax.coastlines()
# ax.legend()
ax.set_global()
ax.gridlines(crs=ccrs.PlateCarree(), linewidth=0.5, color='black',
    draw_labels=True, alpha=0.5, linestyle='--')
c = ax.contourf(lon, lat, D_grid, bounds_c, norm=col.LogNorm(), transform=data_crs, extend = 'both')

plt.colorbar(c, cax=cbar_ax, ax=ax)


plt.show()



# plt.figure()
# plt.title('Diffusivity')
# s = plt.contourf(D_grid, cmap = cm.PuBu_r )
# plt.colorbar(s)
#
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
# # # fit normal and laplace distributions
# x = np.arange(-180,181,1)
# # # a_fit, b_fit = sc.laplace.fit(np.asarray(L_angles) * 180/np.pi)
# norm_a, norm_b = sc.norm.fit(np.asarray(L_ang)) # see if we can make this better
# # # y = sc.laplace(scale = b_fit)
# y_norm = sc.norm(norm_a, norm_b)
#
# # #%%
#
# plt.figure()
# plt.ylabel('frequency')
# plt.xlabel(r'Angle ($\theta$)')
# plt.hist(np.asarray(L_ang), 360, density = True, stacked =True) # check if normalization is done correctly
# plt.plot(x, y_norm.pdf(x), color = 'r')
# plt.show()
#
# # plt.plot(x, y.pdf(x), color = 'r')
#
# # plt.legend(['Fitted normal distribution', 'Data'])
# #
# # plt.show()
# #
# s, p = sc.normaltest(L_ang)
#
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