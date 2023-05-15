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
    def __init__(self, fname_Cell, fname_MPS, fname_Thetas):
        
        self.fname_Cell     = fname_Cell
        self.fname_MPS      = fname_MPS
        self.fname_Thetas   = fname_Thetas

        self.Thetas_loaded = False
        

    def load_Cell(self):
        self.figname_Cell       = self.fname_Cell[self.fname_Cell.index("_")+1:].strip(".txt")
        self.ell, self.C_ell    = np.loadtxt(DATA_PATH + self.fname_Cell, unpack=True)


    def load_matter_PS(self):
        self.figname_matterPS   = self.fname_MPS[self.fname_MPS.index("_")+1:].strip(".txt")
        self.k_h_Mpc, self.Pk   = np.loadtxt(DATA_PATH + self.fname_MPS, skiprows=1, unpack=True)
        self.k_eq               = np.loadtxt(DATA_PATH + self.fname_MPS, max_rows=1) * u.h / u.Mpc 

        self.k_h_Mpc *= u.h / u.Mpc




    def load_Thetas(self):
        self.figname_Thetas       = self.fname_Thetas[self.fname_Thetas.index("_")+1:].strip(".txt")
        self.Thetas_loaded = True 
        file = DATA_PATH+self.fname_Thetas
        Nells   = len(np.loadtxt(file, skiprows=2, max_rows=1))
        cols    = np.arange(1, Nells)

        self.ell_list = np.loadtxt(file, skiprows=1, max_rows=1, usecols=np.arange(2,Nells+1), dtype=int)


        self.eta0   = np.loadtxt(file, skiprows=0, usecols=1, max_rows=1) * u.m
        k           = np.loadtxt(file, skiprows=2, usecols=0) / u.m
        self.k_eta0 = k*self.eta0

        self.Thetas = np.loadtxt(file, skiprows=2, usecols=cols, unpack=True)

        self.Theta_squared_over_k = (np.abs(self.Thetas)**2 / k).to(u.Mpc)
         


    def plot_Cell(self):
        self.load_Cell()
        ylabel = r"$\ell(\ell+1)C_\ell/2\pi$"
        figname = "C_ell_" + self.figname_Cell
        plot.plot_C_ell(self.ell, self.C_ell,
                        fname=figname,
                        ylabel=ylabel, ypad=10,
                        save=SAVE, temp=TEMP, push=PUSH)
        
    def comp_Cell(self, file_compare):
        ell_compare, C_ell_compare = np.loadtxt(DATA_PATH + file_compare, unpack=True) 
        ylabel = r"$\ell(\ell+1)C_\ell/2\pi$"
        plot.comp_C_ell(self.ell, self.C_ell - C_ell_compare,
                        ylabel=ylabel, ypad=10,
                        save=SAVE, temp=TEMP, push=PUSH)

    def plot_matter_power_spectrum(self):
        self.load_matter_PS()
        ylabel = r"$P(k)\quad[(\mathrm{Mpc/h})^3]$"
        xlabel = r"$k \quad [\mathrm{h/Mpc}]$"
        fname = "Matter_PS_" + self.figname_matterPS
        plot.plot_matter_PS(self.k_h_Mpc, self.Pk, self.k_eq,
                            fname=fname,
                            ylabel=ylabel, 
                            xlabel=xlabel, 
                            ypad=10,
                            save=SAVE, temp=TEMP, push=PUSH)


    def plot_Thetas(self, ells=np.array([])):
        if not self.Thetas_loaded:
            self.load_Thetas() 

        res = np.in1d(ells, self.ell_list)
        if not res.all():
            invalid = ells[res==False] 
            print("Attempting to plot ell-values not contained in the dataset.")
            print(f"  ell={invalid} not in dataset.")
            exit()
        else:

            ell_idx = np.where(np.isin(self.ell_list, ells))[0] 
            Thetas = self.Thetas[ell_idx]
            ell_string = "_l"
            for l in ells:
                ell_string += "_" + str(l)

        figname = "Theta_" + self.figname_Thetas + ell_string
        ylabel=r"$\Theta_\ell(k)$"
        plot.plot_Thetas(self.k_eta0, Thetas, ells, fname=figname,
                         ylabel=ylabel, logy=False, xlim=[5e-1, 5e2],
                         save=SAVE, temp=TEMP, push=PUSH)
        

    def plot_Integrand(self, ells=np.array([])):
        if not self.Thetas_loaded:
            self.load_Thetas() 

        res = np.in1d(ells, self.ell_list)
        if not res.all():
            invalid = ells[res==False] 
            print("Attempting to plot ell-values not contained in the dataset.")
            print(f"  ell={invalid} not in dataset.")
            exit()
        else:

            ell_idx = np.where(np.isin(self.ell_list, ells))[0] 
            integrands = self.Theta_squared_over_k[ell_idx]
            ell_string = "_l"
            for l in ells:
                ell_string += "_" + str(l)

        figname = "Theta_squared_over_k_" + self.figname_Thetas + ell_string
        ylabel=r"$|\Theta_\ell(k)|^2 / k \quad[\mathrm{Mpc}]$"
        plot.plot_Thetas(self.k_eta0, integrands, ells, fname=figname,
                         ylabel=ylabel, logy=False, ypad=10, xlim=[5e-1, 1e2],
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
f1_Cell="cells_nx5000_nk5443_nlogk10886.txt" 
f1_MPS="matterPS_nk1000.txt"
f1_Thetas="thetas_nk2000_nx5000.txt"

pspec = PowerSpectrum(f1_Cell, f1_MPS, f1_Thetas)

# SAVE=True
# PUSH=True
# TEMP=True

pspec.plot_Cell()
pspec.plot_matter_power_spectrum()
pspec.plot_Thetas(ells=np.array([2, 4, 7, 10]))
pspec.plot_Integrand(ells=np.array([2, 4, 7, 10]))


