# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 10:23:54 2020

@author: Alex
"""

import numpy as np
import matplotlib.pyplot as plt
#Import the NLL_1 function which  calculates the NLL without the cross section term
from Analysis_Methods import NLL_1

x0 = 0.5
x1 = 0.65
x2 = 0.73

#Define a routine which calculated the curvature at the minimum of a lagrange polynomial interpolated
#through 3 sample points
def d2dx(x0,x1, x2,x3, f):
    y0, y1, y2 = f(x0), f(x1), f(x2)
    t0 = (y0*2)/((x0-x1)*(x0-x2))
    t1 = (y1*2)/((x1-x0)*(x1-x2))
    t2 = (y2*2)/((x2-x0)*(x2-x1))
    return (t0+t1+t2)
    
#Defin a routine which in 1D interpolates a lagrange polynomial through 3 samples, and/
#then calculated the minima of the interpolated parabola
def par_min(x0, x1, x2, f):
    y0, y1, y2 = f(x0), f(x1), f(x2)
    #convergence depends on a residue calculation, so a control bit is used
    #so that the residue isnt referenced before the variable is declared
    control = 0
    #initialise residue as something higher than the threshold convergence value
    res = 10
    while res > 0.00000001:
        #Algorithm will produce zero-division if the sample y values are identical, so make sure
        #loop ends before this is computed
        y0, y1, y2 = f(x0), f(x1), f(x2)
        if y0==y2 or y0==y1 or y1==y2:
            break
        else:
            #x3 is the minima for the interpolated lagrange poly
            x3 = 0.5*((x2**2-x1**2)*y0+(x0**2-x2**2)*y1+(x1**2-x0**2)*y2)/\
                ((x2-x1)*y0+(x0-x2)*y1+(x1-x0)*y2)
            a = np.array([y0,y1,y2,f(x3)])
            b = np.argmax(a)
            #figure out which sample point is the highest, and remove this from the sample set!
            if control == 1:
                #calculate residue to end loop when convergence complete
                res = np.sqrt((x3-xm)**2)
            c = [x0, x1, x2, x3]
            df = d2dx(c[0], c[1], c[2], c[3], f)
            del c[b]
            c = np.array(c)
            c.sort()
        
            x0, x1, x2 = c
            #print(res)
            control = 1
            xm = x3
    std_dev = 1/np.sqrt(df)
    return x3, std_dev

#Definte a function that turns the NLL into a one dimensional function, for ease of plotting
#Parabolic minimiser here is coded as a function of one variable
def f2(x):
    return NLL_1([x, 2.4e-3])

#Define some sample points near the minima, guessed by eye by plotting the function
x0, x1, x2 = 0.6, 0.7, 0.75


print("the function minimum and standard deviation are", par_min(x0, x1, x2, f2))

#Sample the parameter space for plotting
pi = np.linspace(0, np.pi/4, 100)
u_v = []
for i in pi:
    u_v.append([i, 2.4*1e-3])
    
#Compute the NLL for the sampled theta values
NLL_v = []
for i in u_v:
    NLL_v.append(NLL_1(i))
    
#Define pyplot parameters so plots hopefully render well on any python build        
plt.figure(figsize=(10,6))
plt.rc('xtick', labelsize=15)   
plt.rc('ytick', labelsize=15)    
plt.rc('figure', titlesize=45)
plt.rc('axes', titlesize=20)     # fontsize of the axes title
plt.rc('axes', labelsize=20)
plt.rc('legend', fontsize=22)

#Plot NLL
plt.plot(pi, NLL_v)
plt.grid()

plt.xlabel('θ (radians)')
plt.ylabel('NLL')
plt.legend(['NLL(θ)']) 
plt.show()

#Using the computed minima, define functions that scan upwards until the NLL increases
#and decreases by 0.5, correpsonding to one standard deviation
theta_min = 0.6763585471865108
NLL_min = f2(theta_min)

def scan_up():
    diff = 0
    dx = 0
    while diff < 0.5:
        new_NLL = f2(theta_min+dx)
        diff = new_NLL-NLL_min
        dx += 0.0001
    return dx

def scan_down():
    diff = 0
    dx = 0
    while diff < 0.5:
        new_NLL = f2(theta_min+dx)
        diff = new_NLL-NLL_min
        dx -= 0.0001
    return dx


print("Standard deviation in positive direction is", scan_up())
print("Standard deviation in negative direction is", scan_down())
    
    

    

