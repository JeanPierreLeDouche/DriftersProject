# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 09:44:38 2020

@author: Gebruiker
"""
import numpy as np
import cmath

point1 = np.asarray([1,1])
point2 = np.asarray([1,-1])
point3 = np.asarray([-1,-1])
point4 = np.asarray([-1,1])
point0 = np.asarray([0,0])
point5 = np.asarray([-1,0])
point6 = np.asarray([1,0])
point7 = np.asarray([0,-1])
point8 = np.asarray([1,0])

vec1 = point7 - point0
vec2 = point1 - point0

def angle_between_vecs(vector1, vector2):
    x1 = vector1[0]
    y1 = vector1[1]
    
    x2 = vector2[0]
    y2 = vector2[1]    
        
    dot = x1*x2 + y1*y2      # dot product between [x1, y1] and [x2, y2]
    det = x1*y2 - y1*x2      # determinant
    angle = np.arctan2(det, dot) * 180 / np.pi  # angle in degrees [-180, 180]
    
    return angle

angle = angle_between_vecs(vec1, vec2)
print(angle)