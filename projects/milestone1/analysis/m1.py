import numpy as np 
import matplotlib.pyplot as plt 
import astropy.units as u 
from astropy.constants import c  
import os 
import warnings 
warnings.filterwarnings("ignore", category=DeprecationWarning )

"""
All calculations and such are performed in this file,
but the actual plotting is done in plot.py
"""
import plot

DATA_PATH = "/home/vetle/Documents/master_studies/subjects/V23/AST5220/projects/milestone1/data/"


global SAVE
global TEMP
global PUSH
SAVE = False # If False, the figures produced are only displayed, not saved. 
TEMP = False  
PUSH = False

class ExpansionHistory:

    def __init__(self, 
                 fname_cosmology,
                 time_unit=u.Myr,
                 length_unit=u.Mpc,
                 SN_burn=1000):
        self.data   = self.load(fname_cosmology) 
        self.x      = self.data[0]
        
        self.time_unit = time_unit
        self.length_unit = length_unit
        
        self.load_important_times(fname_cosmology) 
        self.load_H_parameters()

        self.load_eta_and_t()
        self.load_omegas() 

        self.load_supernova_data()
        self.load_supernovafit(SN_burn)

    def load(self, filename, skiprows=1):
        data = np.loadtxt(DATA_PATH + filename, 
                          unpack=True, 
                          skiprows=skiprows)
        return data 
    
    def load_supernova_data(self, filename="supernovadata.txt"):
        self.supernovadata = self.load(filename)
        self.N_SN_data = len(self.supernovadata[0])

        
    def load_important_times(self, filename):
        important_times = np.loadtxt(DATA_PATH + filename, 
                                    unpack=True, 
                                    max_rows=1)
        
        equality_times_x = important_times[0:3]
        equality_times_t = important_times[3:6]
        equality_times_z = np.exp(-equality_times_x) - 1
        self.x_mr_eq, self.x_ml_eq, self.x_acc_onset = equality_times_x 
        self.z_mr_eq, self.z_ml_eq, self.z_acc_onset = equality_times_z
        
        t_mr_eq, t_ml_eq, t_acc_onset = equality_times_t
        self.t_mr_eq = (t_mr_eq*u.s).to(u.yr)
        self.t_ml_eq = (t_ml_eq*u.s).to(u.Gyr)
        self.t_acc_onset = (t_acc_onset*u.s).to(u.Gyr)

        t0 = important_times[6] * u.s 
        eta0 = important_times[7] * u.m 
        self.t0 = t0.to(u.Gyr)
        self.eta0_over_c = (eta0/c).to(u.Gyr)


    

    def load_H_parameters(self, convert_units=True):
        # Load Hp, Hp'(x) and H''(x)
        H_params = self.data[1:4] / u.s 

        if convert_units: 
            H_params = H_params.to(100*u.km / u.s / u.Mpc) 

        self.Hp, self.dHp_dx, self.ddHp_ddx = H_params 


    def load_eta_and_t(self, convert_units=True):
        self.eta_t_loaded = True 
        self.eta = self.data[4] * u.m 
        self.t   = self.data[5] * u.s 

        if convert_units:
            self.eta = self.eta.to(self.length_unit)
            self.t = self.t.to(self.time_unit)
        else:
            pass 

    def load_omegas(self):
        self.omegas_loaded = True 
        OmegaB, OmegaCDM, OmegaLambda = self.data[6:9]
        OmegaR, OmegaNu, OmegaK = self.data[9:12]
    
        self.OmegaM = OmegaB + OmegaCDM
        self.OmegaRel = OmegaR + OmegaNu 
        self.OmegaLambda = OmegaLambda



    def luminosity_distance(self, data):
        ### Return z and dL for a given data file
        #  limited to region with observational data 
        z_sim    =  data[0]
        dL_sim   = (data[-1]*u.m).to(u.Gpc)

        z_data = self.supernovadata[0]
        z_min = np.min(z_data) - 0.005
        z_max = np.max(z_data) + 0.01

        # z_sim = self.z_of_x(x_sim) 
        z_sim_mask = (z_sim <= z_max) & (z_sim >= z_min)
        z_sim = z_sim[z_sim_mask]

        dL_sim = dL_sim[z_sim_mask]
        
        return z_sim, dL_sim 




    def z_of_x(self, x):
        return np.exp(-x) - 1    
    

    def plot_eta_H_over_c(self):
        ### Plot eta*Hp/c ### 
    
        ylabel = r'$\eta \mathcal{H} / c $'
        eta_Hp_over_c = (self.eta*self.Hp / c).to(1)

        plot.plot_single_param(self.x, eta_Hp_over_c, 
                               "compare_eta_H_over_c.pdf", 
                               self.x_mr_eq, self.x_ml_eq, 
                               acc=self.x_acc_onset,
                               legend=True,
                               xlabel=r"x", ylabel=ylabel, 
                               xlim=[-16,2], ylim=[0.75, 4], 
                               yticks=[1,2,3,4],
                               log=False, 
                               save=SAVE, temp=TEMP, push=PUSH)

    def plot_dH_ddH_over_H(self):
        ### Compare H'/H and H''/H with analytical approx. ### 

        dHp_over_Hp = self.dHp_dx / self.Hp
        ddHp_over_Hp = self.ddHp_ddx / self.Hp

        dH_label = r"$\frac{\mathcal{H}'(x)}{\mathcal{H}(x)}$"
        ddH_label = r"$\frac{\mathcal{H}''(x)}{\mathcal{H}(x)}$"


        plot.compare_dH_and_ddH_over_H(self.x, dHp_over_Hp, ddHp_over_Hp, 
                                       dH_label, ddH_label, 
                                       self.x_mr_eq, self.x_ml_eq, 
                                       title=None, 
                                       save=SAVE, temp=TEMP, push=PUSH)
    
    
    def plot_Hp(self):
        ### Plot Hp(x) ### 

        ylabel = r'$\mathcal{H}\: \left[ \frac{100\,\mathrm{km/s}}{\mathrm{Mpc}} \right] $'
        
        plot.plot_single_param(self.x, self.Hp, "compare_Hp.pdf", 
                                self.x_mr_eq, self.x_ml_eq, self.x_acc_onset,
                                xlabel=r"$x$", ylabel=ylabel, 
                                xlim=[-16,3], ylim=[1e-1, 1e4], 
                                legend=True, legendloc='lower left', 
                                save=SAVE, temp=TEMP, push=PUSH)
    

    def plot_eta_t(self):
        ### Plot eta and t ### 

        ylabel = rf"Time $[\mathrm{{{self.time_unit}}}]$"
        eta_over_c = (self.eta / c).to(self.time_unit)
        plot.plot_t_and_eta(self.x, self.t, eta_over_c, 
                            fname="t_and_eta_c.pdf", 
                            mr_eq=self.x_mr_eq, 
                            mL_eq=self.x_ml_eq,
                            ylabel=ylabel, ylim=[1e-11, 5e6],
                            acc=self.x_acc_onset, 
                            save=SAVE, temp=TEMP, push=PUSH)

    def plot_omegas(self):
        ### Plot density parameters ### 
        title = r"$\Omega_i(x)$"

        plot.plot_omega_params(self.x, 
                               self.OmegaM, 
                               self.OmegaRel, 
                               self.OmegaLambda, 
                               fname="omega_i_of_x.pdf", 
                               mr_eq=self.x_mr_eq, 
                               ml_eq=self.x_ml_eq,
                               acc=self.x_acc_onset, 
                               xlabel=r'$x$',
                               title=title, 
                               save=SAVE, temp=TEMP, push=PUSH)


    def plot_dL(self, fname_dL_planck, fname_dL_fitted):
        ### Plot dL from data and compare with simulation with Planck parameters  
        ### plot_fit=True: Plot dL from the simulation with parameters obtained from supernovafit.
        #       Planck results included for comparison
        
        z, dL, dL_error = self.supernovadata
        planck_data     = self.load(fname_dL_planck)
        fit_data        = self.load(fname_dL_fitted)

        z_planck, dL_planck = self.luminosity_distance(planck_data)
        z_fit, dL_fit       = self.luminosity_distance(fit_data)

        data    = [z,        dL, dL_error]
        planck  = [z_planck, dL_planck.value]
        fit     = [z_fit,    dL_fit.value]


        plot.plot_dL(data, planck, fit, fname="dL_z_compare_fitted.pdf", save=SAVE, temp=TEMP)




    def load_supernovafit(self, burn):
        # Load result from supernova fit 
        # Leave out the first N=burn results  
        supernova_mcmc_results = self.load("supernovafit.txt", skiprows=1+burn)
        self.chi2, self.h_SN, self.OmegaM_SN, self.OmegaK_SN = supernova_mcmc_results

        self.OmegaLambda_SN = 1 - self.OmegaM_SN - self.OmegaK_SN
        self.H0_SN = self.h_SN * (100 * u.km / u.s / u.Mpc)
        
        self.bestfit_idx = np.argmin(self.chi2)
        self.chi2_min = np.min(self.chi2) 
        self.chi2_1sigma = self.chi2 < self.chi2_min + 3.53
        self.chi2_2sigma = self.chi2 < self.chi2_min + 8.02 




        self.best_fit = list(supernova_mcmc_results[:,self.bestfit_idx])
        self.best_fit[0] /= self.N_SN_data
        self.best_fit[1] *=  100 #* u.km / u.s / u.Mpc



    def plot_supernova_fit_omegas(self, burn=1000):
        # Plot OmegaM-OmegaK confidence region.
        # chi2, h, OmegaM, OmegaK = self.load_supernovafit(burn)
    
        # OmegaLambda = 1 - OmegaM - OmegaK 
        # chi2min = np.min(chi2)
        # chi2_1sigma = chi2 < chi2min + 3.53
        # chi2_2sigma = chi2 < chi2min + 8.02

        plot.plot_OmegaM_OmegaLambda_plane(self.OmegaM_SN, self.OmegaLambda_SN, 
                                        self.chi2_1sigma, self.chi2_2sigma, chi2_min=self.bestfit_idx,
                                        fname=f"mcmc_supernova_fit_Nburn{burn}.pdf", 
                                        save=SAVE, temp=TEMP)

    def plot_supernova_fit_H0_pdf(self, burn=1000):
        # Plot H0 pdf
        # h = self.load_supernovafit(burn)[1]
        # bestfit_idx = np.argmin(self.load_supernovafit)

        # H0 = h * (100*u.km / u.s / u.Mpc)

        H0_mean = np.mean(self.H0_SN)
        H0_std  =  np.std(self.H0_SN) 
        H0_var  =  np.var(self.H0_SN) 

        
        H0_min = H0_mean - 4*H0_std 
        H0_max = H0_mean + 4*H0_std 
        N_bins = 100 
        bins   = np.linspace(H0_min, H0_max, N_bins)

        H0_gaussian_distr = (2*np.pi*H0_std**2)**(-1/2) * np.exp(-(bins - H0_mean)**2 / (2 * H0_var))

        plot.plot_H0_posterior_pdf(self.H0_SN, bins, H0_gaussian_distr,
                                    fname=f'H0_pdf_Nburn{burn}.pdf', 
                                    bestfit = self.best_fit[1],
                                    save=SAVE, temp=TEMP)

    def supernova_bestfit_table(self):
        # Plot H0 pdf
        # chi2, h, OmegaM, OmegaK = self.load_supernovafit(burn)
        # H0 = h * (100*u.km / u.s / u.Mpc)

        # best_fit_idx = np.argmin(chi2)
        # chi2_best = chi2[best_fit_idx]

        
        # chi2best_over_N = chi2_best/self.N_SN_data
        # best_fit = [H0[best_fit_idx].value, OmegaM[best_fit_idx], OmegaK[best_fit_idx], chi2best_over_N]

        standard_deviations = [np.std(self.H0_SN.value), np.std(self.OmegaM_SN), np.std(self.OmegaK_SN), None]
        plot.bestfit_supernova_table(self.best_fit, standard_deviations, save=SAVE)

    def supernova_plot_chi2(self, burn=1000):
        chi2_over_N = self.load_supernovafit(burn)[0] / self.N_SN_data
        chi2_mean = np.mean(chi2_over_N)
        chi2_std  = np.std(chi2_over_N) 
        chi2_var  = np.var(chi2_over_N) 

        
        chi2_min = chi2_mean - 4*chi2_std 
        chi2_max = chi2_mean + 4*chi2_std 
        N_bins = 100 
        bins   = np.linspace(chi2_min, chi2_max, N_bins)

        plot.plot_fit_posterior_pdf(chi2_over_N, bins, mean=chi2_mean,
                                    fname="chi2_pdf.pdf", 
                                    param_name=r"$\chi^2 / N$", xlabel=r"$\chi^2 / N$",
                                    save=SAVE, temp=TEMP)



    def best_fit_values(self, burn=1000):
        # chi2, h, OmegaM, OmegaK = self.load_supernovafit(burn)
        # H0 = h * (100*u.km / u.s / u.Mpc)

        # bestfit = np.argmin(np.abs(chi2/self.N_SN_data - 1))
        print(np.min(self.H0_SN));exit()
        H0_min = self.H0_SN[self.chi2_1sigma].min()
        H0_max = self.H0_SN[self.chi2_1sigma].max() 
        print(self.H0_SN[self.bestfit_idx] - H0_min)
        print(self.H0_SN[self.bestfit_idx] - H0_max)
        print(np.std(self.H0_SN))

        # print(OmegaM[bestfit])
        # print(OmegaK[bestfit])

        
        exit()
        H0_std = np.std(H0)

        OmegaM_std = np.std(OmegaM)
        OmegaK_std = np.std(OmegaK)


    def make_table(self):
        """
        Funker ikke. Pandas er no mig innimellom...
        """

        # print(self.t_mr_eq)
        # print(self.t_acc_onset)
        # print(self.t_ml_eq)
        mr_eq = [self.x_mr_eq, self.z_mr_eq, self.t_mr_eq]
        ml_eq = [self.x_ml_eq, self.z_ml_eq, self.t_ml_eq]
        acc_onset = [self.x_acc_onset, self.z_acc_onset, self.t_acc_onset]
        # print(self.t0,"\n", self.eta0)
        # exit()
        plot.time_table(mr_eq, ml_eq, acc_onset, self.t0, self.eta0_over_c,
                        save=SAVE, temp=TEMP, push=PUSH)

Cosmology = ExpansionHistory("cosmology_new.txt")

SAVE = True 
# TEMP = True 

### COMPLETE ### 
# Cosmology.plot_dH_ddH_over_H()
# Cosmology.plot_eta_H_over_c()
# Cosmology.plot_Hp()
# Cosmology.plot_eta_t()
# Cosmology.plot_omegas()
# Cosmology.plot_dL(fname_dL_fitted="bestfit_dL.txt", fname_dL_planck="planck_dL.txt")
# Cosmology.plot_supernova_fit_omegas()
# Cosmology.plot_supernova_fit_H0_pdf()

Cosmology.supernova_bestfit_table()
### CHECK ### 


# Cosmology.best_fit_values()
# Cosmology.supernova_plot_chi2()

