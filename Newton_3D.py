# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 10:23:54 2020

@author: Alex
"""

import numpy as np
import matplotlib.pyplot as plt

#Import relevant methods from central library
from Analysis_Methods import NLL_2, data1, data2, data3, L




    


#Calculate the 3 dimensional partial derivatives using the chain rule, firstly 
#we calculate the rate that λ changes with respect to each parameter
def d2λ_dθ2(u, i):
    a = (1.267*L)/data3[i]
    return -8*data3[i]*u[2]*data2[i]*(np.sin(a*u[1])**2)*np.cos(4*u[0])

def d2λ_dm2(u, i):
    a = (1.267*L)/data3[i]
    return -2*data3[i]*u[2]*(a**2)*data2[i]*(np.sin(2*u[0])**2)*np.cos(2*a*u[1])

def d2λ_dk2(u, i):
    return 0

def d2λ_dθdm(u, i):
    a = (1.267*L)/data3[i]
    return -2*data3[i]*u[2]*a*data2[i]*np.sin(2*a*u[1])*np.sin(4*u[0])

def d2λ_dkdθ(u, i):
    a = (1.267*L)/data3[i]
    return -2*data3[i]*data2[i]*(np.sin(a*u[1])**2)*np.sin(4*u[0])

def d2λ_dkdm(u, i):
    a = (1.267*L)/data3[i]
    return -data3[i]*data2[i]*(np.sin(2*u[0])**2)*a*np.sin(2*a*u[1])




def dλ_dθ(u, i):
    a = (1.267*L)/data3[i]
    return -2*data3[i]*u[2]*data2[i]*(np.sin(a*u[1])**2)*np.sin(4*u[0])

def dλ_dm(u, i):
    a = (1.267*L)/data3[i]
    return -data3[i]*u[2]*data2[i]*(np.sin(2*u[0])**2)*a*np.sin(2*a*u[1])

def dλ_dk(u, i):
    a = (1.267*L)/data3[i]
    return data3[i]*data2[i]*(1-(np.sin(a*u[1])**2)*(np.sin(2*u[0])**2))


#Using the chain rule, we calculate the change in NLL by summing over the rate of change with
#respect to each λ, and multilpying by the partial derivs of λ w.r.t each parameter
def d2f_dθ2(u):
    sum = 0
    for i in range(0, 200):
        if data1[i] == 0:
            sum += d2λ_dθ2(u, i)
        elif data1[i] != 0:
            a = (1.267*L)/data3[i]
            b = np.sin(a*u[1])**2
            c = np.sin(2*u[0])**2
            λ = u[2]*data3[i]*data2[i]*(1-b*c)
            
            x1 = d2λ_dθ2(u, i)
            x2 = (1-(data1[i]/λ))
            x3 = dλ_dθ(u, i)
            x4 = (data1[i]/(λ**2))
            sum += (x1*x2) + (x3**2)*x4
    return sum

def d2f_dm2(u):
    sum = 0
    for i in range(0, 200):
        if data1[i] == 0:
            sum += d2λ_dm2(u, i)
        elif data1[i] != 0:
            a = (1.267*L)/data3[i]
            b = np.sin(a*u[1])**2
            c = np.sin(2*u[0])**2
            λ = u[2]*data3[i]*data2[i]*(1-b*c)
            
            x1 = d2λ_dm2(u, i)
            x2 = (1-(data1[i]/λ))
            x3 = dλ_dm(u, i)
            x4 = (data1[i]/(λ**2))
            sum += (x1*x2) + (x3**2)*x4
    return sum

def d2f_dk2(u):
    sum = 0
    for i in range(0, 200):
        if data1[i] == 0:
            sum += d2λ_dk2(u, i)
        elif data1[i] != 0:
            a = (1.267*L)/data3[i]
            b = np.sin(a*u[1])**2
            c = np.sin(2*u[0])**2
            λ = u[2]*data3[i]*data2[i]*(1-b*c)
            
            x1 = d2λ_dθ2(u, i)
            x2 = (1-(data1[i]/λ))
            x3 = dλ_dk(u, i)
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
            λ = u[2]*data3[i]*data2[i]*(1-b*c)
            
            x1 = d2λ_dθdm(u, i)
            x2 = (1-(data1[i]/λ))
            x3 = dλ_dθ(u, i)
            x4 = dλ_dm(u, i)
            x5 = (data1[i]/(λ**2))
            sum += (x1*x2) + (x3*x4)*x5
    return sum

def d2f_dkdθ(u):
    sum = 0
    for i in range(0, 200):
        if data1[i] == 0:
            sum += d2λ_dkdθ(u, i)
        elif data1[i] != 0:
            a = (1.267*L)/data3[i]
            b = np.sin(a*u[1])**2
            c = np.sin(2*u[0])**2
            λ = u[2]*data3[i]*data2[i]*(1-b*c)
            
            x1 = d2λ_dkdθ(u, i)
            x2 = (1-(data1[i]/λ))
            x3 = dλ_dk(u, i)
            x4 = dλ_dθ(u, i)
            x5 = (data1[i]/(λ**2))
            sum += (x1*x2) + (x3*x4)*x5
    return sum

def d2f_dkdm(u):
    sum = 0
    for i in range(0, 200):
        if data1[i] == 0:
            sum += d2λ_dkdm(u, i)
        elif data1[i] != 0:
            a = (1.267*L)/data3[i]
            b = np.sin(a*u[1])**2
            c = np.sin(2*u[0])**2
            λ = u[2]*data3[i]*data2[i]*(1-b*c)
            
            x1 = d2λ_dkdm(u, i)
            x2 = (1-(data1[i]/λ))
            x3 = dλ_dk(u, i)
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
            λ = u[2]*data3[i]*data2[i]*(1-b*c)
            
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
            λ = u[2]*data3[i]*data2[i]*(1-b*c)
            
            x1 = dλ_dm(u, i)
            x2 = (1-(data1[i]/λ))
            
            sum += (x1*x2)
    return sum

def df_dk(u):
    sum = 0
    for i in range(0, 200):
        if data1[i] == 0:
            sum += dλ_dk(u, i)
        elif data1[i] != 0:
            a = (1.267*L)/data3[i]
            b = np.sin(a*u[1])**2
            c = np.sin(2*u[0])**2
            λ = u[2]*data3[i]*data2[i]*(1-b*c)
            
            x1 = dλ_dk(u, i)
            x2 = (1-(data1[i]/λ))
            
            sum += (x1*x2)
    return sum




#Compute hessian matrix!
def hess(u):
    return np.array([[d2f_dθ2(u), d2f_dθdm(u), d2f_dkdθ(u)], [d2f_dθdm(u), d2f_dm2(u), d2f_dkdm(u)], [d2f_dkdθ(u), d2f_dkdm(u), d2f_dk2(u)]])

#Compute gradient vector!
def grad(u):
    return np.array([df_dθ(u), df_dm(u), df_dk(u)])

#Calculate newton method and convergence criteria!
def Newton(u_in):
    u_1 = u_in
    res = 10
    for i in range(10):
        g = grad(u_1)
        h = hess(u_1)
        inv = np.linalg.inv(h)
        u_2 = u_1 - np.matmul(inv, g)
        res = np.linalg.norm(u_2-u_1)
        u_1 = u_2
    print(u_2)
    return u_2


#Start the initial guess at the previous values for theta and mass, with a guess of the cross section
#parameter alpha

u_in = [0.73607156, 0.00211372, 3]



#Compute the newton method and return the converged minima
u = Newton(u_in)
