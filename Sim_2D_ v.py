# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 10:23:54 2020

@author: Alex
"""

import numpy as np
import matplotlib.pyplot as plt

#To validate the newton 2D method, hardcode an example function
def f(u):
    return -np.exp(-u[0]**3/3 + u[0] -u[1]**2)

#Hardcode the analytic partial derivatives so that the grad vector
# and hessian matrix can be calculated
def d2f_dy2(u):
    return f(u)*(((-2*u[1])**2)-2)

def d2f_dx2(u):
   return f(u)*(((-u[0]**2+1)**2)-2*u[0])

def d2f_dxdy(u):
    return f(u)*(-u[0]**2+1)*(-2*u[1])

def df_dx(u):
    return f(u)*(-u[0]**2+1)

def df_dy(u):
    return f(u)*(-2*u[1])

#Calculate the hessian matrix/curvature matrix
def hess(u):
    return np.array([[d2f_dx2(u), d2f_dxdy(u)], [d2f_dxdy(u), d2f_dy2(u)]])

#Calculate the gradient vector
def grad(u):
    return np.array([df_dx(u), df_dy(u)])

#Impliment the newton method algorithm
def Newton(u_in):
    u_1 = u_in
    #use a residual method to break the loop once we have converged to a value
    res = 10
    while res > 0:
        g = grad(u_1)
        h = hess(u_1)
        inv = np.linalg.inv(h)
        u_2 = u_1 - np.matmul(inv, g)
        res = np.linalg.norm(u_2-u_1)
        u_1 = u_2
    return u_2

#Define an initial guess
u_in = (0.6, 0.0021172686496694995)

#Compute newton algo
u = Newton(u_in)

print(" The curvature matrix is ", hess(u))
print(" The function minimum is ", u)

#VISUALISE 2D PARABOLA LIKE EXPONENTIAL FUNCTION

#Import some extra data visualisation libraries

from pylab import *
from mpl_toolkits.mplot3d import Axes3D

#define meshgrid that we will now plot
x = np.arange(-1.5, 2.5, 0.3)
y = np.arange(-1.5, 2.5, 0.3)
X, Y = np.meshgrid(x, y)
Z = -np.exp(-X**3/3 + X -Y**2)

#chose a nice colour map, here lighter colours corresponding to lower function values
cmap = 'gnuplot2_r'
fig = figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z, cmap= cmap)
title(cmap)