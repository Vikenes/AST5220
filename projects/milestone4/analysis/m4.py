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
M3_PATH = "/home/vetle/Documents/master_studies/subjects/V23/AST5220/projects/milestone3/data/"
DATA_PATH = "/home/vetle/Documents/master_studies/subjects/V23/AST5220/projects/milestone4/data/"
TEST_PATH = DATA_PATH + "testruns/"
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
    def __init__(self, f1_Cell, f1_MPS):
        
        # f1_Cell=f1_Cell[f1_Cell.index("_")+1:]
        # print(f1_Cell)
        self.f1_Cell = f1_Cell[f1_Cell.index("_")+1:].strip(".txt")
        self.f1_MPS = f1_MPS[f1_MPS.index("_")+1:].strip(".txt")
        self.ell, self.C_ell = np.loadtxt(DATA_PATH + f1_Cell, unpack=True)
        self.k, self.Pk      = np.loadtxt(DATA_PATH + f1_MPS , unpack=True) 
        # self.ell, self.C_ell1 = np.loadtxt(DATA_PATH + f1, unpack=True)
        # self.ell2, self.C_ell2 = np.loadtxt(TEST_PATH + f2, unpack=True)
        # self.ell3, self.C_ell3 = np.loadtxt(TEST_PATH + f3, unpack=True)

        # self.plot()
        # self.plot_diff(self.C_ell2, self.C_ell3)

    def plot_Cell(self):
        ylabel = r"$\ell(\ell+1)C_\ell/2\pi$"
        fname = "C_ell_" + self.f1_Cell
        plot.plot_C_ell(self.ell, self.C_ell,
                        fname=fname,
                        ylabel=ylabel, ypad=10,
                        save=SAVE, temp=TEMP, push=PUSH)
        
    def comp_Cell(self, file_compare):
        ell_compare, C_ell_compare = np.loadtxt(DATA_PATH + file_compare, unpack=True) 
        ylabel = r"$\ell(\ell+1)C_\ell/2\pi$"
        plot.comp_C_ell(self.ell, self.C_ell - C_ell_compare,
                        ylabel=ylabel, ypad=10,
                        save=SAVE, temp=TEMP, push=PUSH)

    def plot_matter_power_spectrum(self):
        ylabel = r"$P(k)\quad[(\mathrm{Mpc/h})^3]$"
        xlabel = r"$k \quad [\mathrm{h/Mpc}]$"
        fname = "Matter_PS_" + self.f1_MPS
        plot.plot_matter_PS(self.k, self.Pk,
                            fname=fname,
                            ylabel=ylabel, 
                            xlabel=xlabel, 
                            ypad=10,
                            save=SAVE, temp=TEMP, push=PUSH)


    def plot_source_function(self, xlim=XLIM, ylim=None, no=0, fname=None):
        # plot.K_FIGNAMES = self.k_fnames
        # if hasattr(self, 'source'):
            # pass
        # else:
            # self.load_source_func()
        
        if no == 0:
            if fname is None:
                fname = "S_alone"
            plot.plot_source_function(self.x, self.S_tilde,
                                      k_legends=self.k_labels,
                                      ylabel=r"$\tilde{S}$",
                                      fname=fname,
                                      xlim=xlim,ylim=ylim,
                                      save=SAVE, push=PUSH, temp=TEMP)

        if no == 1:
            if fname is None:
                fname = "S_j5"
            plot.plot_source_function(self.x, self.S_tilde_j5, 
                                      k_legends=self.k_labels,
                                      ylabel=r"$\tilde{S} j_5 $",
                                      fname=fname,
                                      xlim=xlim,ylim=ylim,
                                      save=SAVE, push=PUSH, temp=TEMP)

        elif no == 2:
            if fname is None:
                fname = "S_j50"
            plot.plot_source_function(self.x, self.S_tilde_j50, 
                                      k_legends=self.k_labels,
                                      ylabel=r"$\tilde{S} j_{50} $",
                                      fname=fname,
                                      xlim=xlim,ylim=ylim,
                                      save=SAVE, push=PUSH, temp=TEMP)
            # plt.title(r"$\tilde{S}(k,x) j_{50}$")
            # plt.plot(self.x, self.S_tilde_j50[2])
            # plt.xlim(xlim)
            # plt.ylim(ylim)
            # plt.show()
        elif no == 3:
            if fname is None:
                fname = "S_j500"
            plot.plot_source_function(self.x, self.S_tilde_j500, 
                                      k_legends=self.k_labels,
                                      ylabel=r"$\tilde{S} j_{500} $",
                                      fname=fname,
                                      xlim=xlim,ylim=ylim,
                                      save=SAVE, push=PUSH, temp=TEMP)
            # plt.title(r"$\tilde{S}(k,x) j_{500}$")
            # plt.plot(self.x, self.S_tilde_j500[2])
            # plt.xlim(xlim)
            # plt.ylim(ylim)
            # plt.show()
        else:
            pass     

   


# pspec = PowerSpectrum(f1="thisworks.txt", f2="cellspar.txt", f3="cellsnewx.txt")
# pspec = PowerSpectrum(f1_Cell="cells.txt", f1_MPS="matterPS.txt")#, f3="cellsignoreterms.txt")
pspec = PowerSpectrum(f1_Cell="cells_nx5000_nk5443_nlogk10886.txt", f1_MPS="matterPS_nk1000.txt")

SAVE=True
PUSH=True
pspec.plot_Cell()
TEMP=True
pspec.plot_Cell()

# pspec.comp_Cell(file_compare="thisworks.txt")
# pspec.plot_matter_power_spectrum()

