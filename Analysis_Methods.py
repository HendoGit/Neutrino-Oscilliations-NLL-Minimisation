# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 10:00:02 2020

@author: Alex
"""
import numpy as np

data = np.loadtxt('data.txt')

n = np.int(len(data)/2)

data1 = data[0:n]
data2 = data[n:2*n]
data3 = []

for i in range(0, 200):
    data3.append(0.025+i*0.05)
data3 = np.array(data3)


#s: θ23 =π/4 , ∆m =2.4 × 10−3   L = 295.
L = 295

def prob(E, u):
    a = np.sin(2*u[0])**2
    b = np.sin((1.267*u[1]*L)/E)**2
    return 1-a*b

    
def λ_1(u):
    Pvals = prob(data3,u)
    data4 = []
    for i in range(0,200):
        data4.append(data2[i]*Pvals[i])
    return data4
    
    
    
def NLL_1(u):
    sum = 0
    λ_ = λ_1(u)
    for i in range(0, 200):
        if λ_[i] != 0 and data1[i] == 0:
            sum += λ_[i]
        else:
            LL = λ_[i] + data1[i]*(np.log(data1[i]/λ_[i]) - 1) 
            sum += LL
        
    return sum


def λ2(u):
    Pvals = prob(data3,u)
    data4 = []
    for i in range(0,200):
        data4.append(u[2]*data3[i]*data2[i]*Pvals[i])
    return data4
    
    
#Input to NLL in the form [θ, mass]
def NLL_2(u):
    sum = 0
    λ_ = λ2(u)
    for i in range(0, 200):
        if λ_[i] != 0 and data1[i] == 0:
            sum += λ_[i]
        else:
            LL = λ_[i] + data1[i]*(np.log(data1[i]/λ_[i]) - 1) 
            sum += LL
        
    return sum