# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 13:17:12 2020

@author: Gebruiker
"""

import numpy as np

def k_p2_func(k_xx, k_xy, k_yy):
    
    theta = 0.5 * np.arctan( 2 * k_xy / (k_xx - k_yy)) # eq 9 Rühs
    k_p2 = k_xx * np.sin(theta)**2 - k_xy * np.sin(2*theta) + k_yy * np.cos(theta)**2 # eq 8 Rühs
       
    return k_p2

def k_davis(v_res, d_res):
    # v_res, d_res are scalar values of the residual velocity and displacement
    k = -1 * v_res * d_res 
    
    return k

def k_disp(d, d_2, prev_d, prev_d_2, delta_t): # displacement arguments all residual
    prev_s = prev_d * prev_d_2# eq 7 Rühs
    s = d * d_2    # eq 7  
    
    delta_s = s - prev_s  
    k = 0.5 * delta_s / delta_t # eq 6 Rühs 
    
    return k

def K(k_davis_p2, k_disp_p2):
    K = (k_davis_p2 + k_disp_p2)/2
    return K

v_x = 1. # being residual velocity in x
v_y = 1. # res. vel. in y

d_y = 1. # res. disp. in y
d_x = 1. # res. disp in x

d_y_old = 1. # prev res. disp. in y
d_x_old = 1. # prev res. disp. in x

K_davis = k_p2_func(k_davis(v_x, d_x), k_davis(v_x, d_y), k_davis(v_y, d_y))
K_disp = k_p2_func(k_disp(d_x, d_x, d_x_old, d_x_old, 24*3600 ), k_disp(d_x, d_y, d_x_old, d_y_old, 24*3600 ), k_disp(d_y, d_y, d_y_old, d_y_old, 24*3600, 24*3600  )
K = K(K_davis, K_disp)
                   