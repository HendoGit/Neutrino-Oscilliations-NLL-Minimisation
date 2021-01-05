# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 10:23:54 2020

@author: Alex
"""

import numpy as np
import matplotlib.pyplot as plt
from Analysis_Methods import NLL_1
from Univariate import univariate

#Import univariate method and apply it to mininise mass first
#Extract the function minima, the curvature, and the standard deviation estimate
extract = univariate(0.5, 0.6,0.7 ,0.0021, 0.0023, 0.0025, NLL_1, c=1)

#Print the final result for this method
print("theta, del_m = ", extract[0:2])
print("sig_theta, sig_m = ", extract[4:])

        







