# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 10:23:54 2020

@author: Alex
"""

import numpy as np
import matplotlib.pyplot as plt

#Import NLL routine for 2 input params as well as input data needed for hardcoded derivs
from Analysis_Methods import NLL_1, data1, data2, data3, L 




#Compute the analytical partial derivatives so that the gradient and hessian can be
#calculated
def d2λ_dm2(u, i):
    a = (1.267*L)/data3[i]
    return -2*(a**2)*data2[i]*(np.sin(2*u[0])**2)*np.cos(2*a*u[1])

def d2λ_dθ2(u, i):
    a = (1.267*L)/data3[i]
    return -8*data2[i]*(np.sin(a*u[1])**2)*np.cos(4*u[0])

def d2λ_dθdm(u, i):
    a = (1.267*L)/data3[i]
    return -2*a*data2[i]*np.sin(2*a*u[1])*np.sin(4*u[0])

def dλ_dm(u, i):
    a = (1.267*L)/data3[i]
    return -data2[i]*(np.sin(2*u[0])**2)*a*np.sin(2*a*u[1])

def dλ_dθ(u, i):
    a = (1.267*L)/data3[i]
    return -2*data2[i]*(np.sin(a*u[1])**2)*np.sin(4*u[0])




#Using the chain rule, the derivative of the NLL with respect to the params is split
#into the rate of change of λ with respect to the parameter, and the rate of change of 
#NLL with respect to λ
#NOTE we then sum over all of the λ to get the analtic partial derivatives!!!
def d2f_dm2(u):
    sum = 0
    for i in range(0, 200):
        if data1[i] == 0:
            sum += d2λ_dm2(u, i)
        elif data1[i] != 0:
            a = (1.267*L)/data3[i]
            b = np.sin(a*u[1])**2
            c = np.sin(2*u[0])**2
            λ = data2[i]*(1-b*c)
            
            x1 = d2λ_dm2(u, i)
            x2 = (1-(data1[i]/λ))
            x3 = dλ_dm(u, i)
            x4 = (data1[i]/(λ**2))
            sum += (x1*x2) + (x3**2)*x4
    return sum

def d2f_dθ2(u):
    sum = 0
    for i in range(0, 200):
        if data1[i] == 0:
            sum += d2λ_dθ2(u, i)
        elif data1[i] != 0:
            a = (1.267*L)/data3[i]
            b = np.sin(a*u[1])**2
            c = np.sin(2*u[0])**2
            λ = data2[i]*(1-b*c)
            
            x1 = d2λ_dθ2(u, i)
            x2 = (1-(data1[i]/λ))
            x3 = dλ_dθ(u, i)
            x4 = (data1[i]/(λ**2))
            sum += (x1*x2) + (x3**2)*x4
    return sum

def d2f_dθdm(u):
    sum = 0
    for i in range(0, 200):
        if data1[i] == 0:
            sum += d2λ_dθdm(u, i)
        elif data1[i] != 0:
            a = (1.267*L)/data3[i]
            b = np.sin(a*u[1])**2
            c = np.sin(2*u[0])**2
            λ = data2[i]*(1-b*c)
            
            x1 = d2λ_dθdm(u, i)
            x2 = (1-(data1[i]/λ))
            x3 = dλ_dθ(u, i)
            x4 = dλ_dm(u, i)
            x5 = (data1[i]/(λ**2))
            sum += (x1*x2) + (x3*x4)*x5
    return sum

def df_dθ(u):
    sum = 0
    for i in range(0, 200):
        if data1[i] == 0:
            sum += dλ_dθ(u, i)
        elif data1[i] != 0:
            a = (1.267*L)/data3[i]
            b = np.sin(a*u[1])**2
            c = np.sin(2*u[0])**2
            λ = data2[i]*(1-b*c)
            
            x1 = dλ_dθ(u, i)
            x2 = (1-(data1[i]/λ))
            
            sum += (x1*x2)
    return sum

def df_dm(u):
    sum = 0
    for i in range(0, 200):
        if data1[i] == 0:
            sum += dλ_dm(u, i)
        elif data1[i] != 0:
            a = (1.267*L)/data3[i]
            b = np.sin(a*u[1])**2
            c = np.sin(2*u[0])**2
            λ = data2[i]*(1-b*c)
            
            x1 = dλ_dm(u, i)
            x2 = (1-(data1[i]/λ))
            
            sum += (x1*x2)
    return sum



#Define hessian matrix!
def hess(u):
    return np.array([[d2f_dθ2(u), d2f_dθdm(u)], [d2f_dθdm(u), d2f_dm2(u)]])

#Define gradient vector!
def grad(u):
    return np.array([df_dθ(u), df_dm(u)])

#Define newton method, with residual calculated each iteration so that loop breaks when converged!
def Newton(u_in):
    u_1 = u_in
    res = 10
    while res > 0:
        g = grad(u_1)
        h = hess(u_1)
        inv = np.linalg.inv(h)
        u_2 = u_1 - np.matmul(inv, g)
        res = np.linalg.norm(u_2-u_1)
        u_1 = u_2
    return u_2

#Start at previous guess from parabolic methods!
u_in = (0.7346504620818721, 0.0021172686496694995)


u = Newton(u_in)

cov = np.linalg.inv(hess(u))
err_theta = np.sqrt(cov[0,0])
err_m = np.sqrt(cov[1,1])
print("The estimate for the fit parameters = ", u)
print("The estimate for the fit errors = ", err_theta, err_m)
