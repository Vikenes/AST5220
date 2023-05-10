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

DATA_PATH = "/home/vetle/Documents/master_studies/subjects/V23/AST5220/projects/milestone4/data/"

global SAVE 
global PUSH
global TEMP
global XLIM
global MREQ
global TCEND 
SAVE        = False 
PUSH        = False
TEMP        = False 
XLIM        = [-14.5, 0]
MREQ        = True 
TCEND       = False 

class PowerSpectrum:
    def __init__(self, f1, f2, f3):
        


        self.ell, self.C_ell = np.loadtxt(DATA_PATH + f1, unpack=True)
        self.ell2, self.C_ell2 = np.loadtxt(DATA_PATH + f2, unpack=True)
        self.ell3, self.C_ell3 = np.loadtxt(DATA_PATH + f3, unpack=True)

        plt.plot(self.ell,  self.C_ell, c='blue')
        plt.plot(self.ell2, self.C_ell2, '--', c='k')
        plt.plot(self.ell3, self.C_ell3, ':', alpha=0.5)
        plt.xscale('log')
        # plt.yscale('log')
        plt.show()



pspec = PowerSpectrum(f1="thisworks.txt", f2="cells.txt", f3="cellspar.txt")

# SAVE=True
# PUSH=True
# TEMP=True

