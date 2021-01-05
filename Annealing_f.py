# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 10:23:54 2020

@author: Alex
"""

import numpy as np
import matplotlib.pyplot as plt

#import the simulated annealing methods and the 3 dimensional NLL function
from Annealing import Main
from Analysis_Methods import NLL_2
    
#apply the sim annealing method and calculate the parameter chain
chain = Main(NLL_2)

#caculate the steps for plotting
steps = np.arange(len(chain))

#format graphs
plt.figure(figsize=(10,6))
plt.rc('xtick', labelsize=15)   
plt.rc('ytick', labelsize=15)    
plt.rc('figure', titlesize=45)
plt.rc('axes', titlesize=20)     # fontsize of the axes title
plt.rc('axes', labelsize=20)
plt.rc('legend', fontsize=22)

#plot the variation of each parameter as a function of step
plt.plot(steps, ParamChain[:, 0], 'r-')
plt.plot(steps, ParamChain[:, 1], 'b-')
plt.plot(steps, ParamChain[:, 2], 'g-')

plt.grid()
plt.xlabel("No. of Steps") 
plt.ylabel("Parameter Value") 
plt.legend(["θ","Δm", "α"]) 
plt.show()

#plot the variation of system temperature as a function of step
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
array([0.77812787, 0.00217278, 1.82031335])



