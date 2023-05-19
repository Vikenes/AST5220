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
# EXTDATA_PATH = "/home/vetle/Documents/master_studies/subjects/V23/AST5220/projects/milestone4/data/external/"
COMPARISON_PATH = DATA_PATH + "comparison_data/"
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
        C_ell_data              = np.loadtxt(DATA_PATH + self.fname_Cell, unpack=True)
        self.ell                = C_ell_data[0]
        self.C_ell              = C_ell_data[1]
        self.C_ell_components   = C_ell_data[2:] 

    def load_planck_C_ell(self, fname):
        CMB_planck = np.loadtxt(COMPARISON_PATH + fname, unpack=True)
        ell_planck = CMB_planck[0]
        C_ell_planck = CMB_planck[1]
        error_planck = CMB_planck[2:4]
        return ell_planck, C_ell_planck, error_planck
    
    def load_Pk_external(self, fname):
        data = np.loadtxt(COMPARISON_PATH + fname, unpack=True)
        return data 


    def load_matter_PS(self):
        self.k_h_Mpc, self.Pk   = np.loadtxt(DATA_PATH + self.fname_MPS, skiprows=1, unpack=True)
        self.k_eq               = np.loadtxt(DATA_PATH + self.fname_MPS, max_rows=1) * u.h / u.Mpc 

        self.k_h_Mpc *= u.h / u.Mpc




    def load_Thetas(self):
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
         


    def plot_Cell(self, fname_planck):
        ell_planck, C_ell_planck, error_planck = self.load_planck_C_ell(fname_planck)

        self.load_Cell()
        ylabel = r"$\ell(\ell+1)C_\ell/2\pi$"
        figname = self.fname_Cell
        plot.plot_C_ell(self.ell, self.C_ell,
                        ell_planck, C_ell_planck,
                        error_planck,
                        fname=figname,
                        ylabel=ylabel, ypad=10,
                        logy=False, fill=True,
                        save=SAVE, temp=TEMP, push=PUSH)
        

    def comp_Cell(self, file_compare):
        ell_compare, C_ell_compare = np.loadtxt(DATA_PATH + file_compare, unpack=True) 
        ylabel = r"$\ell(\ell+1)C_\ell/2\pi$"
        plot.comp_C_ell(self.ell, self.C_ell - C_ell_compare,
                        ylabel=ylabel, ypad=10,
                        save=SAVE, temp=TEMP, push=PUSH)

    def plot_Cell_components(self, fname_planck):
        ell_planck, C_ell_planck, error_planck = self.load_planck_C_ell(fname_planck)

        self.load_Cell()
        ylabel = r"$\ell(\ell+1)C_\ell/2\pi$"
        figname = self.fname_Cell.strip(".txt")
        plot.plot_C_ell_components(self.ell, self.C_ell, self.C_ell_components,
                        ell_planck, C_ell_planck,
                        error_planck,
                        fname=figname,
                        ylabel=ylabel, ypad=10,
                        logy=False, fill=True,
                        save=SAVE, temp=TEMP, push=PUSH)
        




    def plot_matter_power_spectrum(self, fname_galaxy_survey, fname_wmap):
        galaxy_data = self.load_Pk_external(fname_galaxy_survey)
        wmap_data = self.load_Pk_external(fname_wmap)
        wmap_data[2] = np.abs(wmap_data[1] - wmap_data[2])

        self.load_matter_PS()

        ylabel = r"$P(k)\quad[(\mathrm{Mpc/h})^3]$"
        xlabel = r"$k \quad [\mathrm{h/Mpc}]$"
        fname = self.fname_MPS.strip(".txt")
        plot.plot_matter_PS(self.k_h_Mpc, self.Pk, self.k_eq,
                            galaxy_data,
                            wmap_data,
                            fname=fname,
                            ylabel=ylabel, 
                            xlabel=xlabel, 
                            xlim=[3e-4, 5e-1], ylim=[2e2, 5e4],
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

        figname = self.fname_Thetas.split(".txt")[0] + ell_string
        ylabel=r"$\Theta_\ell(k)$"


        plot.plot_Thetas(self.k_eta0, Thetas, ells, fname=figname,
                         ylabel=ylabel, logy=False, logx=False, xlim=[5e-1, 5e2],
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

        figname = "integrand_" + self.fname_Thetas.split(".txt")[0] + ell_string
        ylabel=r"$|\Theta_\ell(k)|^2 / k \quad[\mathrm{Mpc}]$"
        plot.plot_Thetas(self.k_eta0, integrands, ells, fname=figname,
                         ylabel=ylabel, logy=False, logx=False, ypad=10, xlim=[5e-1, 1e2],
                         save=SAVE, temp=TEMP, push=PUSH)


   

fname_Cell="cells_components_nx1500_nk10886_nlogk10886.txt"
fname_MPS="matterPS_nk1000.txt"
fname_Thetas="thetas_nx1500_nk10886_nlogk10886.txt"

fname_planck = "planck_cell_low.txt"
fname_galaxy_survey = "reid_DR7.txt"
fname_wmap = "wmap_act.txt"

pspec = PowerSpectrum(fname_Cell, fname_MPS, fname_Thetas)

SAVE=True
# PUSH=True
# TEMP=True

pspec.plot_Cell_components(fname_planck)
pspec.plot_matter_power_spectrum(fname_galaxy_survey, fname_wmap)
pspec.plot_Thetas(ells=np.array([6, 100, 200]))
pspec.plot_Integrand(ells=np.array([2, 4, 7, 10]))


