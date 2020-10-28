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
import matplotlib.colors as col

def k_p2_func(k_xx, k_xy, k_yy):
    theta = 0.5 * np.arctan(2 * k_xy / (k_xx - k_yy))  # eq 9 Rühs
    k_p2 = k_xx * np.sin(theta) ** 2 - k_xy * np.sin(2 * theta) + k_yy * np.cos(theta) ** 2  # eq 8 Rühs

    return k_p2


def k_davis(v_res, d_res):
    # v_res, d_res are scalar values of the residual velocity and displacement
    k = -1 * v_res * d_res

    return k


def k_disp(d, d_2, prev_d, prev_d_2, delta_t):  # displacement arguments all residual
    prev_s = prev_d * prev_d_2  # eq 7 Rühs
    s = d * d_2  # eq 7

    delta_s = s - prev_s
    k = 0.5 * delta_s / delta_t  # eq 6 Rühs

    return k


def K(k_davis_p2, k_disp_p2):
    K = (k_davis_p2 + k_disp_p2) / 2
    return K

data = pickle.load(
    open(r'C:\Users\Ruben\Documents\CLPH\MAIO\ruhsdata.p', "rb"))

r_e = 6.37e6 #m

# verschillende k grids

lat_points = 180
lon_points = 360
K_grid_dav = [[[] for x in range(lon_points)] for x in range(lat_points + 1)]
K_grid_disp = [[[] for x in range(lon_points)] for x in range(lat_points + 1)]
K_grid = [[[] for x in range(lon_points)] for x in range(lat_points + 1)]

buoy_IDs = np.unique(data["ID"])
print(f'There are {len(buoy_IDs)} unique buoys')


#averages to calculate residual velocities in m/s

ve_avg = np.mean(np.asarray(data["VE(CM/S)"]))/100
vn_avg = np.mean(np.asarray(data["VN(CM/S)"]))/100

print(ve_avg, vn_avg)
#per buoy berekeningen
###### Loop over boeien en selecteer de goede__________________________________________________________________
for j, ID in enumerate(buoy_IDs[:]):
    print('Buoy ',j)
    buoy = data.loc[data["ID"] == buoy_IDs[j]]


    if len(buoy['Lon']) < 81:
        print(f'Buoy {j} does not have enough data')
        continue

    if np.asarray(buoy["VE(CM/S)"])[0] >400:
        print(f'Initial speed incorrect, skipping to next buoy')
        continue

    if j in [508, 512]:
        continue

    buoy = buoy[:-1]
    print(f'The buoy length is: {len(buoy)}')
#########

    #Davis method___________________________________________
    # set t_0
    t_0 = len(buoy['Lat'] - 1)

    #Load data
    lons = np.asarray(buoy["Lon"])
    lats = np.asarray(buoy["Lat"])
    ve = np.asarray(buoy["VE(CM/S)"])/100
    vn = np.asarray(buoy["VN(CM/S)"])/100

    #Initialize vectors
    x_displacement = np.zeros(20)
    y_displacement = np.zeros(20)
    x_displacement_avg = np.zeros(20)
    y_displacement_avg = np.zeros(20)
    x_exp = np.zeros(20)
    y_exp = np.zeros(20)
    x_disp_res = np.zeros(20)
    y_disp_res = np.zeros(20)

    i = 0
    v_x_res = ve[t_0 - 1] - ve_avg
    v_y_res = vn[t_0 - 1] - vn_avg

    for t in range(t_0, t_0 - 80, -4):

        x_exp[i] = (i+1) * (3600 * 24) * ve_avg
        x_displacement[i] = np.abs((lons[t_0 - 1] - lons[t - 5]) * np.cos(lats[t-5]*(np.pi/180)) * (r_e/360))
        x_disp_res[i] = x_displacement[i] - x_exp[i]

        y_exp[i] = (i + 1) * (3600 * 24) * vn_avg
        y_displacement[i] = np.abs((lats[t_0 - 1] - lats[t - 5]) * (r_e / 360))
        y_disp_res[i] = y_displacement[i] - y_exp[i]

        i += 1

    k_xx_dav = k_davis(v_x_res, x_disp_res)
    k_xy_dav = k_davis(v_x_res, y_disp_res)
    k_yx_dav = k_davis(v_y_res, x_disp_res)
    k_yy_dav = k_davis(v_y_res, y_disp_res)

    k_xy_dav = (k_xy_dav + k_yy_dav)/2

    k_p2_dav = k_p2_func(k_xx_dav, k_xy_dav, k_yy_dav)

    # print(k_p2_dav)

    #####Disp methode
    t_0 = 0

    x_displacement = np.zeros(21)
    y_displacement = np.zeros(21)
    x_displacement_avg = np.zeros(21)
    y_displacement_avg = np.zeros(21)
    x_exp = np.zeros(21)
    y_exp = np.zeros(21)
    x_disp_res = np.zeros(21)
    y_disp_res = np.zeros(21)

    i = 0

    dt = 24*3600


    if j == 254:
        print(t_0)

    for t in range(t_0 + 4,t_0+80, 4):

        x_exp[i] = (i + 1) * (3600 * 24) * ve_avg #volgensmij eerder al correct gedefinieerd
        x_displacement[i] = np.abs((lons[t_0] - lons[t]) * np.cos(lats[t] * (np.pi / 180)) * (r_e / 360))
        x_disp_res[i] = x_displacement[i] - x_exp[i]

        y_exp[i] = (i + 1) * (3600 * 24) * vn_avg
        y_displacement[i] = np.abs((lats[t_0 - 1] - lats[t]) * (r_e / 360))
        y_disp_res[i] = y_displacement[i] - y_exp[i]

        i += 1

    k_xx_disp = k_disp(x_disp_res[1:], x_disp_res[1:],x_disp_res[:-1], x_disp_res[:-1],dt)[:-1]
    k_xy_disp = k_disp(x_disp_res[1:], y_disp_res[1:],x_disp_res[:-1], y_disp_res[:-1],dt)[:-1]
    k_yy_disp = k_disp(y_disp_res[1:], y_disp_res[1:],y_disp_res[:-1], y_disp_res[:-1],dt)[:-1]


    k_p2_disp = k_p2_func(k_xx_disp, k_xy_disp, k_yy_disp)

    # print(k_p2_disp)

    ###________________________________________________________

    k_inf_dav = np.abs(np.mean(k_p2_dav[14:]))
    k_inf_disp = np.abs(np.mean(k_p2_disp[14:]))
    k_inf_t = K(k_inf_dav, k_inf_disp)

    # print(k_inf_dav, k_inf_disp, k_inf_t)

    lats_seen = (np.unique(lats))
    for i, _ in enumerate(lats_seen):
        lats_seen[i] = int(lats_seen[i])

    lats_seen = (np.unique(lats_seen))


    lons_seen = (np.unique(lons))
    for i, _ in enumerate(lons_seen):
        lons_seen[i] = int(lons_seen[i])

    lons_seen = (np.unique(lons_seen))



    # print('Seen lats are', lats_seen)
    # print('Seen lons are', lons_seen)



    lon = int(lons[0])
    lat = int(lats[0])

    new_lat = lat + 90
    new_lon = lon

    for lo in lons_seen:
        for la in lats_seen:
            lon = int(lo)
            lat = int(la)
            K_grid[lat][lon].append(k_inf_t)
            K_grid_dav[lat][lon].append(k_inf_dav)
            K_grid_disp[lat][lon].append(k_inf_disp)


k_i = np.zeros((181,360))
k_i_dav = np.zeros((181,360))
k_i_disp = np.zeros((181,360))


for lat in range(0, 181):
    for lon in range(0, 360):
        K_grid_slice = K_grid[lat][lon]
        if len(K_grid_slice) != 0:
            k_i[lat][lon] = np.mean(np.asarray(K_grid[lat][lon]))
            k_i_dav[lat][lon] = np.mean(np.asarray(K_grid_dav[lat][lon]))
            k_i_disp[lat][lon] = np.mean(np.asarray(K_grid_disp[lat][lon]))

        else:
            k_i[lat][lon] = np.NaN
            k_i_dav[lat][lon] = np.NaN
            k_i_disp[lat][lon] = np.NaN

# print((K_grid))


### Plotting_____________________________________________________________________________________________

longit = np.linspace(0,360 , 360)
latit = np.linspace(-90, 90, 181)

data_crs = ccrs.PlateCarree()

fig, ax  = plt.subplots(1,1,figsize=(10,5), subplot_kw={'projection': ccrs.PlateCarree()})

cbar_ax = fig.add_axes([0.9, 0.2, 0.03, 0.7])
fig.subplots_adjust(hspace=0, wspace=0, top=0.90, right=0.8)

bounds_c = [1e2, 5e2, 1e3, 5e3, 1e4, 5e4, 1e5]

ax.set_title(r'K_inf ')
ax.coastlines()
# ax.legend()
ax.set_global()
ax.gridlines(crs=ccrs.PlateCarree(), linewidth=0.5, color='black',
    draw_labels=True, alpha=0.5, linestyle='--')
c = ax.contourf(longit,latit, k_i, bounds_c, norm=col.LogNorm(), transform=data_crs, extend = 'both')

plt.colorbar(c, cax=cbar_ax, ax=ax)

plt.show()