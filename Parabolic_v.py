# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 10:23:54 2020

@author: Alex
"""

import numpy as np
import matplotlib.pyplot as plt

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
        y0, y1, y2 = f(x0), f(x1), f(x2)
        #Algorithm will produce zero-division if the sample y values are identical, so make sure
        #loop ends before this is computed
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
            control = 1
            xm = x3
    return x3, df

#Define a routine which calculated the curvature at the minimum of a lagrange polynomial interpolated
#through 3 sample points
def d2dx(x0,x1, x2,x3, f):
    y0, y1, y2 = f(x0), f(x1), f(x2)
    #this is essentially the second derivative of the above formula
    t0 = (y0*2)/((x0-x1)*(x0-x2))
    t1 = (y1*2)/((x1-x0)*(x1-x2))
    t2 = (y2*2)/((x2-x0)*(x2-x1))
    return (t0+t1+t2)

#define two sample functions for validation of algo
def f(x):
    return 1+x**2

def f2(x):
    return np.exp(x**2)

#define some sample x values, non-symetric so that we dont automatically get 0
x0, x1, x2 = -0.5, 0, 1


print("for function 1, the minimum and curvature are", par_min(x0, x1, x2, f))
print("for function 2, the minimum and curvature are", par_min(x0, x1, x2, f2))


x = np.linspace(-2, 2, 1000)

def f2(x):
    return np.exp(x**2)


        
plt.figure(figsize=(10,6))
plt.rc('xtick', labelsize=15)   
plt.rc('ytick', labelsize=15)    
plt.rc('figure', titlesize=45)
plt.rc('axes', titlesize=20)     # fontsize of the axes title
plt.rc('axes', labelsize=20)
plt.rc('legend', fontsize=22)

plt.plot(x, f(x), 'r-')
plt.plot(x, f2(x), 'b--')
plt.ylim(0, 10)


plt.grid()
#plt.xlabel("x") 
#plt.ylabel("f(x)") 
plt.xlabel('x')
plt.ylabel('y')
plt.legend([r'$f(x) = {x^2} + 1$', r'$f(x) = {e^{x^2}} + 1$']) 
plt.show()



