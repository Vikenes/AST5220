import numpy as np 
import matplotlib.pyplot as plt 
import astropy.units as u 
import os 
import warnings 
warnings.filterwarnings("ignore", category=DeprecationWarning )

from scipy.signal import find_peaks 


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
global XTICKS
SAVE        = False 
PUSH        = False
TEMP        = False 
XLIM        = [-14.5, 0]
XTICKS = [2, 10, 50, 100, 500, 1000, 2000]
class PowerSpectrum:
    def __init__(self, fname_Cell, fname_MPS, fname_Thetas):
        
        self.fname_Cell     = fname_Cell
        self.fname_MPS      = fname_MPS
        self.fname_Thetas   = fname_Thetas

        ### Load CMB power spectrum  
        C_ell_data              = np.loadtxt(DATA_PATH + self.fname_Cell, unpack=True)
        self.ell                = C_ell_data[0]
        self.C_ell              = C_ell_data[1]  
        self.C_ell_components   = C_ell_data[2:] # C_ell from individual terms in source func.  

        ### Load Matter power spectrum 
        self.k_h_Mpc, self.Pk   = np.loadtxt(DATA_PATH + self.fname_MPS, skiprows=1, unpack=True)
        self.k_eq               = np.loadtxt(DATA_PATH + self.fname_MPS, max_rows=1) * u.h / u.Mpc 
        self.k_h_Mpc *= u.h / u.Mpc
        
        ### Load photon multipoles 
        
        # Number of ells used to store Theta_ell
        Nells   = len(np.loadtxt(DATA_PATH + fname_Thetas, skiprows=2, max_rows=1))
        cols    = np.arange(1, Nells)

        # ell's corresponding to theta_ell values  
        self.ell_list = np.loadtxt(DATA_PATH + fname_Thetas, skiprows=1, max_rows=1,usecols=np.arange(2,Nells+1), dtype=int)

        # Horizon today, and the dimless quantity k*eta0
        self.eta0   = np.loadtxt(DATA_PATH + fname_Thetas, skiprows=0, usecols=1, max_rows=1) * u.m
        k           = np.loadtxt(DATA_PATH + fname_Thetas, skiprows=2, usecols=0) / u.m
        self.k_eta0 = k*self.eta0

        # Theta_ell, and the integrand Theta_ell^2 / k
        self.Thetas = np.loadtxt(DATA_PATH + fname_Thetas, skiprows=2, usecols=cols, unpack=True)
        self.Theta_squared_over_k = (np.abs(self.Thetas)**2 / k).to(u.Mpc)





    def load_planck_C_ell(self, fname):
        # Load CMB power spectrum from planck obs 
        CMB_planck = np.loadtxt(COMPARISON_PATH + fname, unpack=True)
        ell_planck = CMB_planck[0]
        C_ell_planck = CMB_planck[1]
        error_planck = CMB_planck[2:4]
        error_planck = np.flip(error_planck, axis=0)
        return ell_planck, C_ell_planck, error_planck
    
    def load_Pk_external(self, fname):
        # Load matter power spectrum from external file 
        data = np.loadtxt(COMPARISON_PATH + fname, unpack=True)
        return data 





    def find_C_ell_extrema(self, find_closest_theta_ell=False, in_units_of_k=False):

        peaks = np.array(find_peaks(self.C_ell)[0])
        troughs = np.array(find_peaks(-self.C_ell)[0][1:])
        ell_peaks = self.ell[peaks]
        ell_troughs = self.ell[troughs]



        if find_closest_theta_ell:
            ell_peaks_closest = self.find_nearest_theta_ell_from_ell_list(ell_peaks)
            ell_troughs_closest = self.find_nearest_theta_ell_from_ell_list(ell_troughs)
            return ell_peaks_closest, ell_troughs_closest
        else:
            if in_units_of_k:
                ell_peaks = (ell_peaks/self.eta0).to(u.Mpc**(-1)) 
                ell_troughs = (ell_troughs/self.eta0).to(u.Mpc**(-1))    
            return ell_peaks, ell_troughs


    def plot_Theta0_at_peaks(self, npeaks=3):
        ell_peaks = self.find_C_ell_extrema()[0][:npeaks]

        data = np.loadtxt(DATA_PATH + "Theta0_of_x_at_peaks.txt", skiprows=1, unpack=True)
        x_dec = np.loadtxt(DATA_PATH + "Theta0_of_x_at_peaks.txt", max_rows=1, usecols=1)
        x = data[0]
        Theta0 = data[1:npeaks+1]
        ylabel = r"$\Theta_0(k,x)$"

        plot.plot_Theta0(x, Theta0, ell_values=ell_peaks,
                            fname="Theta0_at_peaks", x_dec=x_dec,
                            ylabel=ylabel,
                            xlim=[-10, -5], color='red', ypad=10,
                            save=SAVE, push=PUSH, temp=TEMP)


    def plot_Theta0_at_troughs(self, ntroughs=3):
        ell_troughs = self.find_C_ell_extrema()[1][:ntroughs]

        data = np.loadtxt(DATA_PATH + "Theta0_of_x_at_troughs.txt", skiprows=1, unpack=True)
        x_dec = np.loadtxt(DATA_PATH + "Theta0_of_x_at_troughs.txt", max_rows=1, usecols=1)
        x = data[0]
        Theta0 = data[1:ntroughs+1]
        ylabel = r"$\Theta_0(k,x)$"
        
        plot.plot_Theta0(x, Theta0, ell_values=ell_troughs,
                            fname="Theta0_at_troughs", x_dec=x_dec,
                            ylabel=ylabel, color='blue',
                            xlim=[-10, -5],
                            save=SAVE, push=PUSH, temp=TEMP)


    def plot_Theta0_at_peaks_and_troughs(self, npeaks=3, ntroughs=3):
        ell_peaks, ell_troughs = self.find_C_ell_extrema()
        
        x_dec           = np.loadtxt(DATA_PATH + "Theta0_of_x_at_peaks.txt", max_rows=1, usecols=1)

        data_peaks      = np.loadtxt(DATA_PATH + "Theta0_of_x_at_peaks.txt", skiprows=1, unpack=True)
        data_troughs    = np.loadtxt(DATA_PATH + "Theta0_of_x_at_troughs.txt", skiprows=1, unpack=True)
        x               = data_peaks[0]
        Theta0_peaks    = data_peaks[1:npeaks+1]
        Theta0_troughs  = data_troughs[1:ntroughs+1]

        Theta0          = np.concatenate((Theta0_peaks, Theta0_troughs), axis=0)
        # Theta0 -= Theta0[0,0]
        # print(Theta0[:,0]);exit()
        ell             = np.concatenate((ell_peaks[:npeaks], ell_troughs[:ntroughs]))
        ylabel          = r"$\Theta_0(k,x)$"
        
        plot.subplot_Theta0(x, Theta0, ell_values=ell,
                            fname="Theta0_at_peaks_and_troughs", x_dec=x_dec,
                            ylabel=ylabel, ylim=[-0.9, 2.1], legendloc='upper left',
                            xlim=[-11, -5], color='red', ypad=10,
                            save=SAVE, push=PUSH, temp=TEMP)


    def plot_Psi_at_peaks_and_troughs(self, npeaks=3, ntroughs=3):
        ell_peaks, ell_troughs = self.find_C_ell_extrema()
        
        x_dec           = np.loadtxt(DATA_PATH + "Theta0_of_x_at_peaks.txt", max_rows=1, usecols=1)

        data_peaks      = np.loadtxt(DATA_PATH + "Psi_of_x_at_peaks.txt", skiprows=0, unpack=True)
        data_troughs    = np.loadtxt(DATA_PATH + "Psi_of_x_at_troughs.txt", skiprows=0, unpack=True)
        x               = data_peaks[0]
        Psi_peaks    = data_peaks[1:npeaks+1]
        Psi_troughs  = data_troughs[1:ntroughs+1]

        Psi          = np.concatenate((Psi_peaks, Psi_troughs), axis=0)
        ell             = np.concatenate((ell_peaks[:npeaks], ell_troughs[:ntroughs]))
        ylabel          = r"$\Psi(k,x)$"
        
        plot.subplot_Theta0(x, Psi, ell_values=ell,
                            fname="Psi_at_peaks_and_troughs", x_dec=x_dec,
                            ylabel=ylabel, ylim=[-0.8, 0.2], legendloc='upper left',
                            xlim=[-11, -5], color='red', ypad=10, figsize=(15,10),
                            save=SAVE, push=PUSH, temp=TEMP)

    def plot_Psi_plus_Theta0_at_peaks_and_troughs(self, npeaks=3, ntroughs=3):
        ell_peaks, ell_troughs = self.find_C_ell_extrema()
        
        x_dec           = np.loadtxt(DATA_PATH + "Theta0_of_x_at_peaks.txt", max_rows=1, usecols=1)

        Tdata_peaks      = np.loadtxt(DATA_PATH + "Theta0_of_x_at_peaks.txt", skiprows=1, unpack=True)
        Tdata_troughs    = np.loadtxt(DATA_PATH + "Theta0_of_x_at_troughs.txt", skiprows=1, unpack=True)
        x               = Tdata_peaks[0]
        Theta0_peaks    = Tdata_peaks[1:npeaks+1]
        Theta0_troughs  = Tdata_troughs[1:ntroughs+1]

        Theta0          = np.concatenate((Theta0_peaks, Theta0_troughs), axis=0)

        data_peaks      = np.loadtxt(DATA_PATH + "Psi_of_x_at_peaks.txt", skiprows=0, unpack=True)
        data_troughs    = np.loadtxt(DATA_PATH + "Psi_of_x_at_troughs.txt", skiprows=0, unpack=True)
        Psi_peaks    = data_peaks[1:npeaks+1]
        Psi_troughs  = data_troughs[1:ntroughs+1]

        Psi          = np.concatenate((Psi_peaks, Psi_troughs), axis=0)
        Psi_plus_Theta0 = Psi + Theta0

        ell             = np.concatenate((ell_peaks[:npeaks], ell_troughs[:ntroughs]))
        ylabel          = r"$\Psi(k,x)$"
        
        plot.subplot_Theta0(x, Psi_plus_Theta0, ell_values=ell,
                            fname="Psi_at_peaks_and_troughs", x_dec=x_dec,
                            ylabel=ylabel, ylim=[-0.8, 0.9], legendloc='upper left',
                            xlim=[-11, -5], color='red', ypad=10, figsize=(15,10),
                            save=SAVE, push=PUSH, temp=TEMP)
        

    def find_nearest_theta_ell_from_ell_list(self, ell_list):

        ell_closest = []
        for ell in ell_list:
            ell_closest_idx = np.argmin(np.abs(ell - self.ell_list))
            ell_closest.append(self.ell_list[ell_closest_idx])
        return np.array(ell_closest)


    def plot_Cell(self, fname_planck, logx=True, fill=False, n_extrema=0):

        planck = self.load_planck_C_ell(fname_planck)

        ylabel = r"$\ell(\ell+1)C_\ell/2\pi\:(10^6 T_\mathrm{CMB0})^2$"
        figname_split = self.fname_Cell.split("components_")
        figname = figname_split[0] + figname_split[1].strip(".txt")

        plot.plot_C_ell(self.ell, self.C_ell,
                        planck=planck,
                        fname=figname,
                        ylabel=ylabel, ypad=10,
                        xticks=XTICKS, logx=True,
                        logy=False, fill=fill, 
                        save=SAVE, temp=TEMP, push=PUSH)
        

    def plot_Cell_extrema(self, n_extrema=0):

        peaks, troughs = self.find_C_ell_extrema()
        peaks = peaks[:n_extrema]
        troughs = troughs[:n_extrema]

        planck = self.load_planck_C_ell(fname_planck)


        # ylabel = r"$\ell(\ell+1)C_\ell/2\pi\:[\mathrm{\mu K^2}]$"
        ylabel = r"$\ell(\ell+1)C_\ell/2\pi\:\cdot (10^6\, T_\mathrm{CMB0})^2$"
        figname_split = self.fname_Cell.split("components_")
        figname = "peaks_and_troughs_" + figname_split[0] + figname_split[1].strip(".txt")


        plot.plot_C_ell(self.ell, self.C_ell,
                        planck=planck,
                        fname=figname,
                        ylabel=ylabel, ypad=10,
                        logx=True,
                        xticks=XTICKS, 
                        logy=False, fill=False, 
                        peaks=peaks, troughs=troughs,
                        save=SAVE, temp=TEMP, push=PUSH)



    def plot_Cell_components(self):


        ylabel = r"$\ell(\ell+1)C_\ell/2\pi\:\cdot (10^6\, T_\mathrm{CMB0})^2$"

        figname = self.fname_Cell.strip(".txt")
        
        plot.plot_C_ell_components(self.ell, self.C_ell, self.C_ell_components,
                        fname=figname,
                        ylabel=ylabel, ypad=10,
                        logy=False, fill=True, xticks=XTICKS,
                        logx=True,
                        save=SAVE, temp=TEMP, push=PUSH)
        




    def plot_matter_power_spectrum(self, fname_galaxy_survey, fname_wmap):
        fname_MPS_with_Neff = "matterPS_Neff_nk2000.txt"
        _, Pk_Neff   = np.loadtxt(DATA_PATH + fname_MPS_with_Neff, skiprows=1, unpack=True)
        k_eq_Neff = np.loadtxt(DATA_PATH + fname_MPS_with_Neff, max_rows=1) * u.h / u.Mpc 


        galaxy_data = self.load_Pk_external(fname_galaxy_survey)
        wmap_data = self.load_Pk_external(fname_wmap)
        wmap_data[2] = np.abs(wmap_data[1] - wmap_data[2])

        ylabel = r"$P_L(k)\quad[(\mathrm{Mpc/h})^3]$"
        xlabel = r"$k \quad [\mathrm{h/Mpc}]$"
        fname = self.fname_MPS.strip(".txt")
        plot.plot_matter_PS(self.k_h_Mpc, self.Pk, self.k_eq,
                            galaxy_data,
                            wmap_data,
                            Pk_Neff=Pk_Neff, k_eq_Neff=k_eq_Neff,
                            fname=fname,
                            ylabel=ylabel, 
                            xlabel=xlabel, 
                            xlim=[3e-4, 5e-1], ylim=[2e2, 5e4],
                            ypad=10,
                            save=SAVE, temp=TEMP, push=PUSH)


    def plot_Thetas(self, ells=np.array([]), xlim=[5e-1,5e2]):

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
                         ylabel=ylabel, logy=False, logx=False, xlim=xlim,
                         save=SAVE, temp=TEMP, push=PUSH)
        

    def plot_Integrand(self, ells=np.array([]), xlim=[5e-1,5e2]):

        res = np.in1d(ells, self.ell_list)
        if not res.all():
            invalid = ells[res==False] 
            print("Attempting to plot ell-values not contained in the dataset.")
            print(f"  ell={invalid} not in dataset.")
            exit()
        else:

            ell_idx = np.where(np.isin(self.ell_list, ells))[0] 
            normfactor = ells * (ells + 1)
            
            integrands = self.Theta_squared_over_k[ell_idx] * normfactor[:, np.newaxis]
            ell_string = "_l"
            for l in ells:
                ell_string += "_" + str(l)

        figname = "integrand_" + self.fname_Thetas.split(".txt")[0] + ell_string
        ylabel=r"$\ell(\ell+1)|\Theta_\ell(k)|^2 / k \quad[\mathrm{Mpc}]$"
        plot.plot_Thetas(self.k_eta0, integrands, ells, fname=figname,
                         ylabel=ylabel, logy=False, logx=False, ypad=10, xlim=xlim,
                         save=SAVE, temp=TEMP, push=PUSH)


   

# fname_Cell="cells_components_nx1500_nk10886_nlogk10886.txt"
# fname_Thetas="thetas_nx1500_nk10886_nlogk10886.txt"
fname_Cell="cells_components_nx2000_nk21773_nlogk21773.txt"
fname_Thetas="thetas_nx2000_nk21773_nlogk21773.txt"


fname_MPS="matterPS_nk2000.txt"
fname_planck = "planck_cell_low.txt"
fname_galaxy_survey = "reid_DR7.txt"
fname_wmap = "wmap_act.txt"

pspec = PowerSpectrum(fname_Cell, fname_MPS, fname_Thetas)

# SAVE=True
# PUSH=True
# TEMP=True

### CMB power spectrum 
# pspec.plot_Cell(fname_planck, logx=True)
# pspec.plot_Cell_extrema(3)
# pspec.plot_Cell_components()
# pspec.plot_Theta0_at_peaks_and_troughs()

# pspec.plot_Psi_at_peaks_and_troughs()
# pspec.plot_Psi_plus_Theta0_at_peaks_and_troughs()

### Matter power spectrum 
# pspec.plot_matter_power_spectrum(fname_galaxy_survey, fname_wmap)



### Thetas 
pspec.plot_Thetas(ells=[50, 200], xlim=[1e0, 1e3])

### Integrand 
pspec.plot_Integrand(ells=np.array([20, 200, 300]), xlim=[1,5e2])

### NOT USING 
# pspec.plot_Integrand(ells=np.array([2, 5, 10]), xlim=[1, 200])
# pspec.plot_Thetas(ells=[1000, 2000], xlim=[9e2, 2.5e3])

