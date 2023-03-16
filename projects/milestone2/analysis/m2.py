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


data1 = plot.load("recombination.txt", skiprows=2)
# data2 = plot.load("recombination_1e4_1e5.txt", skiprows=2)
# data2 = plot.load("recombination_newode.txt", skiprows=2)

x1 = data1[0]
# xmin = np.abs(x1 + 12).argmin()
# xmax = np.abs(x1).argmin()
# x1 = x1[xmin:xmax]
# data1 = data1[:,xmin:xmax]

# x2 = data2[0]

def load_xe_ne(data):

    Xe = data[1]
    ne = data[2] * u.m**(-3)
    return Xe, ne 

def load_taus(data):
    tau, dtau, ddtau = data[3:6]
    return tau, dtau, ddtau 

def load_gs(data):
    g, dg, ddg = data[6:9]
    return g, dg, ddg 

# g2,dg2,ddg2 = load_gs(data2)

# t2, dt2, ddt2 = load_taus(data2)


def plot_tau():
    t1, dt1, ddt1 = load_taus(data1)
    plt.plot(x1, t1)
    plt.plot(x1, -dt1, '--')
    plt.plot(x1, ddt1, ':')
    plt.xlim(-12,0)
    plt.yscale('log')
    plt.show()

def plot_g():
    g1,dg1,ddg1 = load_gs(data1)
    plt.plot(x1, g1)
    plt.show()
    plt.plot(x1, dg1)#, '--')
    plt.show()
    plt.plot(x1, ddg1)#, ':')
    plt.show()


plot_g()
