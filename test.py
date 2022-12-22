# -*- coding: utf-8 -*-
"""
Created on Thu Dec 22 13:45:42 2022

@author: Lorenzo
"""
import scipy.special as spec
import numpy as np
import scipy.integrate as ints
from math import *
import scipy.optimize as opt
import sys
import matplotlib.pyplot as mpl
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d as plt3d
import matplotlib.lines as ml
import scipy.interpolate as ipt

test1=True
test2=True
file1 = open('Output/MaxLik.dat', 'r')
Lines = file1.readlines()
opt=Lines[-1]
optrue=[1.406 ,  0.330,   0.369,   1.569 ,  0.083 ,  0.010,   1.000 ,   5908.378, 4.0]
tolerance=[0.05,  0.05, 0.05,  0.05 , 0.05 ,  0.05,   0.05 ,  0.05]
mcbest=Lines[-3]
#mcbest=mcbest.split(' ')
resbest=[float(s) for s in mcbest.split()]
resopt=[float(s) for s in opt.split()]
for i in range(0,len(optrue[:-1])):
    val=optrue[i]
    if (abs((val-resopt[i])/val)>=tolerance[i]):
        test1=False
if (abs(resbest[-2]-optrue[-2])>1.0):
        test2=False
        
if(test1):
    print("TEST OPTIMIZATION PASSED")
else:
    print("WARNING: TEST OPTIMIZATION FAILED")
    
if(test2):
    print("TEST BEST FIT PASSED")
else:
    print("WARNING: TEST BEST FIT FAILED")
