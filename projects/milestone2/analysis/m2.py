import numpy as np 
import matplotlib.pyplot as plt 
import astropy.units as u 
import os 
import warnings 
warnings.filterwarnings("ignore", category=DeprecationWarning )

from scipy.integrate import simpson 


"""
All calculations and such are performed in this file,
but the actual plotting is done in plot.py
"""
import plot



save = False # If False, the figures produced are only displayed, not saved. 


data = plot.load("recombination.txt", skiprows=1)
x = data[0]
Xe = data[1]
ne = data[2] * u.m**(-3)
tau, dtau, ddtau = data[3:6]
g, dg, ddg = data[6:9]

gint = simpson(g, x)
# print(gint)
# plt.plot(x, Xe)
# plt.plot(x,   g)
# plt.plot(x,  dg, "--")
# plt.plot(x, ddg, ":")

# plt.plot(x,   tau)
# plt.plot(x, -dtau, "--")
# plt.plot(x, ddtau, ":")
plt.plot(x[ddtau>0], ddtau[ddtau>0])



plt.show()
