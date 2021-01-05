# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 10:23:54 2020

@author: Alex
"""

#Import numerical libraries needed
import numpy as np
import matplotlib.pyplot as plt
import random as rnd

#Set random seeds so that demonstrators can hopefully reproduce identical plots!
iterations = 10000
rnd.seed(0)
np.random.seed(0)

#Implement metropolis hastings algo to choose if we accept or reject new value based on the change in loss
def Metro_hast(delE, T):
    if delE< 0.0:
        return 1
    elif rnd.random() < np.exp(-1.0*delE/T):
        return 1
    return 0

#Implement a proposal function which is bounded -> this feature was added since not bounding 
# the proposal function can result in the function diverging while the temperature is hot!
def prop_fun(u):
    x = u[0] + np.random.normal(0, 0.01)
    if x>5:
        x=5
    if x<-5:
        x=-5
    y = u[1] + np.random.normal(0, 0.01)
    if y>5:
        y=5
    if y<-5:
        y=-5
    z = u[2] + np.random.normal(0, 0.01)
    if z>5:
        z=5
    if z<-5:
        z=-5
    return np.array([x,y,z])


#Define the function that runs the simulated annealing algorithm, and takes a loss function f
#as an input parameter
def Main(f):
    #We define a set of temperatures which decrease exponentially!
    temps = np.arange(0, 15, 0.1)
    temps = 10*np.exp(-temps)
    chain = np.zeros((iterations*len(temps)+1, 3), float, 'C')
    #this eponential temperature decrease empirically from testing gives better convergence!
    chain[0, 0] = 1
    chain[0, 1] = 2
    chain[0, 2] = -1
    
    #we define numpy arrays and a param chain since python lists are slow and this implementation
    #is much faster than any other alternative
    params = np.zeros(3)
    
    #this implementation is adapted from mark scotts, since the concept of using a numpy parameter chain
    #does greatly speed up the computation
    loss = f(chain[0])
    
    for i in range(len(temps)):
        t = temps[i]
        params[0] = chain[iterations*i, 0]
        params[1] = chain[iterations*i, 1]
        params[2] = chain[iterations*i, 2]
        
        loss= f(chain[0])
        
        for j in range(1, iterations+1):
            new_params = prop_fun(params)
            new_loss = f(new_params)
            
            delE = new_loss - loss
            
            res = Metro_hast(delE, t)
            
            if res == 1:
                loss = new_loss
                params = new_params
                
            chain[iterations*i + j, 0] = params[0]
            chain[iterations*i + j, 1] = params[1]
            chain[iterations*i + j, 2] = params[2]
    
    
    return paramChain
    



def f(u):
    return -np.sinc(u[0])*np.sinc(u[1])*np.sinc(u[2])

chain = Main(f)

s = np.arange(len(chain))

plt.figure(figsize=(10,6))
plt.rc('xtick', labelsize=15)   
plt.rc('ytick', labelsize=15)    
plt.rc('figure', titlesize=45)
plt.rc('axes', titlesize=20)     # fontsize of the axes title
plt.rc('axes', labelsize=20)
plt.rc('legend', fontsize=22)


plt.plot(s, chain[:, 0], 'r-')
plt.plot(s, chain[:, 1], 'b-')
plt.plot(s, chain[:, 2], 'g-')

plt.grid()
plt.xlabel("No. of Steps") 
plt.ylabel("Parameter Value") 
plt.legend(["x","y", "z"]) 
plt.show()

temps = np.arange(0, 15, 0.1)
temps = 10*np.exp(-temps)

steps = np.arange(0, 1.5e6, 10000) + 5000

plt.figure(figsize=(10,6))
plt.rc('xtick', labelsize=15)   
plt.rc('ytick', labelsize=15)    
plt.rc('figure', titlesize=45)
plt.rc('axes', titlesize=20)     # fontsize of the axes title
plt.rc('axes', labelsize=20)
plt.rc('legend', fontsize=22)


plt.plot(steps, temps, 'r-')


plt.grid()
plt.xlabel("No. of Steps") 
plt.ylabel("Temperature") 

plt.show()

#VISUALISE 2D MULTIVARIATE SINC FUNCTION

#Import some extra data visualisation libraries

from pylab import *
from mpl_toolkits.mplot3d import Axes3D
x = linspace(-3, 3, 200)
y = linspace(-3, 3, 200)
X,Y = meshgrid(x, y)
Z = -np.sinc(X)*np.sinc(Y)
#chose a nice colour map, here lighter colours corresponding to lower function values
cmap = 'gnuplot2_r'
fig = figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z, cmap= cmap)
title(cmap)




