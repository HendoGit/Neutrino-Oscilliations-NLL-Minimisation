# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 10:23:54 2020

@author: Alex
"""

import numpy as np
import matplotlib.pyplot as plt
from Analysis_Methods import λ2, data, data1, data2, data3




plt.figure(figsize=(10,6))
plt.rc('xtick', labelsize=15)   
plt.rc('ytick', labelsize=15)    
plt.rc('figure', titlesize=45)
plt.rc('axes', titlesize=20)     # fontsize of the axes title
plt.rc('axes', labelsize=20)
plt.rc('legend', fontsize=22)

plt.bar(data3,data1, width=0.1) 


plt.grid()
plt.xlabel("x") 
plt.ylabel("y") 

plt.show()

u = [0.70335073, 0.00218915, 1.42209561]
plt.figure(figsize=(10,6))
plt.rc('xtick', labelsize=15)   
plt.rc('ytick', labelsize=15)    
plt.rc('figure', titlesize=45)
plt.rc('axes', titlesize=20)     # fontsize of the axes title
plt.rc('axes', labelsize=20)
plt.rc('legend', fontsize=22)

plt.bar(data3,λ2(u), width=0.1) 


plt.grid()
plt.xlabel("x") 
plt.ylabel("y") 

plt.show()

u2 = [0.76914108, 0.00217335, 1.85406093]

plt.figure(figsize=(10,6))
plt.rc('xtick', labelsize=15)   
plt.rc('ytick', labelsize=15)    
plt.rc('figure', titlesize=45)
plt.rc('axes', titlesize=20)     # fontsize of the axes title
plt.rc('axes', labelsize=20)
plt.rc('legend', fontsize=22)

plt.bar(data3,λ2(u2), width=0.1) 


plt.grid()
plt.xlabel("x") 
plt.ylabel("y") 

plt.show()

plt.figure(figsize=(10,6))
plt.rc('xtick', labelsize=15)   
plt.rc('ytick', labelsize=15)    
plt.rc('figure', titlesize=45)
plt.rc('axes', titlesize=20)     # fontsize of the axes title
plt.rc('axes', labelsize=20)
plt.rc('legend', fontsize=22)

plt.plot(data3,λ2(u), 'b-') 
plt.plot(data3,λ2(u2), 'r-') 
plt.plot(data3,data1, 'g-') 
plt.xlim(0, 2)
plt.grid()
plt.xlabel("Energy (GeV)") 
plt.ylabel("Number of counts") 
plt.legend(["Parameter Set 1","Parameter Set 2", "Measured data"]) 
plt.show()