# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 10:23:54 2020

@author: Alex
"""

import numpy as np
import matplotlib.pyplot as plt


#Define an adapted parabolic minimiser which operates in one dimension but takes a function of 
#two variables
def par_min(x0, x1, x2, t0, f, c):
    if c == 1:
        y0, y1, y2 = f([x0, t0]), f([x1, t0]), f([x2, t0])
    elif c == 0:
        y0, y1, y2 = f([t0, x0]), f([t0, x1]), f([t0, x2])
    if y0==y2 or y0==y1 or y1==y2:
        raise error
        return
    else:
        
        x3 = 0.5*((x2**2-x1**2)*y0+(x0**2-x2**2)*y1+(x1**2-x0**2)*y2)/((x2-x1)*y0+(x0-x2)*y1+(x1-x0)*y2)
        if c == 1:
            a = np.array([y0,y1,y2,f([x3, t0])])
        elif c == 0:   
            a = np.array([y0,y1,y2,f([t0, x3])])
        b = np.argmax(a)
        c = [x0, x1, x2, x3]
        del c[b]
        x0, x1, x2 = c
        return x3, c

#Define some sample data points around the minima in both the x and y direction
x0, x1, x2 = 0.4, 1.6, 0.9
y0, y1, y2 = -1.1, 0, 1.2

#Define the routine which takes the second derivative of the last parabolic estimate
#by taking the second derivative of the lagrange polynomial
def d2dx(x0,x1, x2, t0, f, c):
    if c == 1:
        y0, y1, y2 = f([x0, t0]), f([x1, t0]), f([x2, t0])
    elif c == 0:
        y0, y1, y2 = f([t0, x0]), f([t0, x1]), f([t0, x2])
    if x0==x2 or x0==x1 or x1==x2:
        raise error
        return
    else:
        t0 = (y0*2)/((x0-x1)*(x0-x2))
        t1 = (y1*2)/((x1-x0)*(x1-x2))
        t2 = (y2*2)/((x2-x0)*(x2-x1))
        return (t0+t1+t2)


#Adapt the previous parabolic algorithm, and now alternate between the first and second 
#variable being the variable that we are minimising
def univariate(x0, x1, x2, y0, y1, y2, f, c):
    #include switch c which controls which variable we minimise first!
    y3 = y1
    x3 = x1
    switch = True
    # c==1 corresponds to minimising x first, then y
    if c == 1:
        #Minimise x first
        x3, x = par_min(x0, x1, x2, y3,f,c=1)
        x0, x1, x2 = x
    #again, use control bit to prevent residue being called before assignment
    control = 0
    res = 10
    #set convergence criteria
    while res > 0.0000000000001:
        #use error handling so that the loop breaks in the event of any errors
        #for example, when near convergence, sample values can become identical, 
        #which will cause a zero-division error, so this prevents that crashing the algo
        try:
            y3, y = par_min(y0, y1, y2, x3,f, c=0)
            y0, y1, y2 = y
            
            dy = d2dx(y0, y1, y2, x3, f, c=0)
        except:
            #break
            pass
        try:
            x3, x = par_min(x0, x1, x2, y3,f, c=1)
            x0, x1, x2 = x
            dx = d2dx(x0, x1, x2, y3, f, c=1)
        except:
            pass
        if control == 1:
             res = np.sqrt((x3-xm)**2 + (y3-ym)**2)
        control = 1
        xm = x3
        ym = y3
    
    sig_x, sig_y = 1/np.sqrt(dx), 1/np.sqrt(dy)
    return x3, y3, dx, dy, sig_x, sig_y


#Defin a sample function to test our algorithm against
def f(u):
    return -np.exp(-u[0]**3/3 + u[0] -u[1]**2)


#Here the data is printed [x_co-orindate, y_co-ordinate, x_curvature, y_curvature]
print("The function minimum and curvature estimates are", univariate(x0, x1, x2, y0, y1, y2, f, c=1)[:4])






