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
                 filename,
                 time_unit=u.Myr,
                 length_unit=u.Mpc):
        self.data   = self.load(filename) 
        self.x      = self.data[0]
        self.x_mr_eq, self.x_ml_eq, self.x_acc_onset = self.load_important_times(filename) 

        self.time_unit = time_unit
        self.length_unit = length_unit

        self.H_loaded = False 
        self.eta_t_loaded = False 
        self.omegas_loaded = False 

    def load(self, filename):
        data = np.loadtxt(DATA_PATH + filename, 
                          unpack=True, 
                          skiprows=1)
        return data 
        
    def load_important_times(self, filename):
        x_mr_eq, x_ml_eq, x_acc_onset = np.loadtxt(DATA_PATH + filename, 
                                                   unpack=True, 
                                                   max_rows=1)
        return x_mr_eq, x_ml_eq, x_acc_onset
    

    def load_H_parameters(self, convert_units=True):
        # Load Hp, Hp'(x) and H''(x)
        self.H_loaded = True 

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

    def z_of_x(self, x):
        return np.exp(-x) - 1    
    

    def plot_eta_H_over_c(self):
        ### Plot eta*Hp/c ### 
        if not self.H_loaded:
            self.load_H_parameters()
        if not self.eta_t_loaded:
            self.load_eta_and_t()
    
        ylabel = r'$\eta \mathcal{H} / c $'

        eta_Hp_over_c = (self.eta*self.Hp / c).to(1)

        plot.plot_single_param(self.x, eta_Hp_over_c, 
                               "compare_eta_H_over_c.pdf", 
                               self.x_mr_eq, self.x_ml_eq, 
                               legend=True,
                               xlabel=r"x", ylabel=ylabel, 
                               xlim=[-16,2], ylim=[0.75, 4], 
                               yticks=[1,2,3,4],
                               log=False, 
                               save=SAVE, temp=TEMP, push=PUSH)

    def plot_dH_ddH_over_H(self):
        ### Compare H'/H and H''/H with analytical approx. ### 
        if not self.H_loaded:
            self.load_H_parameters()


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
        if not self.H_loaded:
            self.load_H_parameters()

        ylabel = r'$\mathcal{H}\: \left[ \frac{100\,\mathrm{km/s}}{\mathrm{Mpc}} \right] $'
        
        plot.plot_single_param(self.x, self.Hp, "compare_Hp.pdf", 
                                self.x_mr_eq, self.x_ml_eq, self.x_acc_onset,
                                xlabel=r"$x$", ylabel=ylabel, 
                                xlim=[-16,3], ylim=[1e-1, 1e4], 
                                legend=True, legendloc='lower left', 
                                save=SAVE, temp=TEMP, push=PUSH)
    

    def plot_eta_t(self):
        ### Plot eta and t ### 
        if not self.eta_t_loaded:
            self.load_eta_and_t()

        # yunit  = f"\mathrm{self.time_unit}" 
        # print(self.time_unit);exit()
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
        if not self.omegas_loaded:
            self.load_omegas()
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

Cosmology = ExpansionHistory("cosmology_new.txt")

# Cosmology.plot_eta_H_over_c()
# Cosmology.plot_dH_ddH_over_H()
# Cosmology.plot_Hp()
# Cosmology.plot_eta_t()
Cosmology.plot_omegas()
exit()


def data(fname="cosmology.txt"):
    """
    Some functions (table()) use data from another file. 
    This is a stupid, but simple way of not having to change all functions for this one purpose. 
    I would have made a class if I had more time.  
    """
    global simulation_data
    global x 
    simulation_data = plot.load(fname)
    x = simulation_data[0]

data()


def load_H_parameters(convert_units=True):
    # Load Hp, Hp'(x) and H''(x)
    H_params = simulation_data[1:4] / u.s 

    if convert_units: 
        H_params = H_params.to(100*u.km / u.s / u.Mpc) 

    Hp, dHp_dx, ddHp_ddx = H_params 

    return Hp, dHp_dx, ddHp_ddx


def load_eta_and_t(convert_units=True):
    eta = simulation_data[4] * u.m 
    t   = simulation_data[5] * u.s 

    if convert_units:
        eta = eta.to(u.Mpc)
        t = t.to(u.Gyr)

    return eta, t 


def load_omegas():
    OmegaB, OmegaCDM, OmegaLambda = simulation_data[6:9]
    OmegaR, OmegaNu, OmegaK = simulation_data[9:12]

    OmegaM = OmegaB + OmegaCDM
    OmegaRel = OmegaR + OmegaNu 

    return OmegaM, OmegaRel, OmegaLambda 


def equality_times(idx=False):
    # Find matter radiation equality and matter-dark energy equality
    # Returns the corresponding x-values by default. 
    # idx=True returns the index of these times. 
    OmegaM, OmegaRel, OmegaLambda = load_omegas()

    Rad_dom_end = np.where(OmegaRel < OmegaM+OmegaLambda)[0][0]
    Lambda_dom_onset = np.where(OmegaLambda > OmegaM+OmegaRel)[0][0]

    MR_eq_idx = np.abs(OmegaRel[:Rad_dom_end+1] - OmegaM[:Rad_dom_end+1]).argmin() 
    ML_eq_idx = np.abs(OmegaM[MR_eq_idx:Lambda_dom_onset+1] - OmegaLambda[MR_eq_idx:Lambda_dom_onset+1]).argmin()
    ML_eq_idx += MR_eq_idx 

        
    if idx:
        return MR_eq_idx, ML_eq_idx
    else:
        return x[MR_eq_idx], x[ML_eq_idx]

def acceleration_onset(idx=False):

    Hp, dHp_dx = load_H_parameters()[0:2]
    a_dot_dot = (np.exp(-x) * Hp * dHp_dx).to((u.Gyr)**(-2)) 

    accleration_onset_idx = np.argmin(np.abs(a_dot_dot))
    if idx:
        return accleration_onset_idx
    else:
        return x[accleration_onset_idx]


def z_of_x(x_):
    return np.exp(-x_) - 1 

def radiation_matter_equality(idx=False):
    OmegaM, OmegaRel, OmegaLambda = load_omegas()
    matter_dom_onset_idx = np.where(OmegaM > OmegaRel+OmegaLambda)[0][0]



def Hp_plot():
    ### Plot Hp(x) ### 
    Hp = load_H_parameters()[0]
    ylabel = r'$\mathcal{H}\: \left[ \frac{100\,\mathrm{km/s}}{\mathrm{Mpc}} \right] $'

    mr_eq, mL_eq = equality_times()
    
    plot.plot_single_param(x, Hp, "compare_Hp.pdf", 
                            mr_eq=mr_eq, mL_eq=mL_eq, acc=acceleration_onset(),
                            xlabel=r"$x$", ylabel=ylabel, 
                            xlim=[-16,3], ylim=[1e-1, 1e4], 
                            legend=True, legendloc='lower left', 
                            save=SAVE, temp=TEMP)



def eta_plot():
    ### Plot eta(x). Not included in report. ### 
    eta = load_eta_and_t()[0]
    ylabel = r'$\eta\:\:[\mathrm{Mpc}]$'
    mr_eq, mL_eq = equality_times()
    plot.plot_single_param(x, eta, "compare_eta.pdf", 
                            xlabel=r"x", ylabel=ylabel, 
                            xlim=[-15,5], ylim=[0.1, 5e5],
                            mr_eq=mr_eq, mL_eq=mL_eq, log=False,
                            legend=True, save=SAVE, temp=TEMP)


def eta_H_plot():
    ### Plot eta*Hp/c ### 
    Hp = load_H_parameters()[0]
    eta = load_eta_and_t()[0]
    ylabel = r'$\eta \mathcal{H} / c $'

    eta_Hp_over_c = (eta*Hp/c).to(1)
    mr_eq, mL_eq = equality_times()
    plot.plot_single_param(x, eta_Hp_over_c, "compare_eta_H_over_c.pdf", 
                            mr_eq=mr_eq, mL_eq=mL_eq, legend=True,
                            xlabel=r"x", ylabel=ylabel, 
                            xlim=[-16,2], ylim=[0.75, 4], 
                            yticks=[1,2,3,4],
                            log=False, save=SAVE, temp=TEMP)


def plot_omegas():
    ### Plot density parameters ### 
    OmegaM, OmegaRel, OmegaLambda = load_omegas()
    title = r"$\Omega_i(x)$"

    plot.plot_omega_params(x, OmegaM, OmegaRel, OmegaLambda, xlabel=r'$x$', 
                        fname="omega_i_of_x.pdf", title=title, save=SAVE, temp=TEMP)
    


def dH_ddH_over_H():
    ### Compare H'/H and H''/H with analytical approx. ### 
    Hp, dHp, ddHp = load_H_parameters()
    m_dom, L_dom = equality_times()

    dH_label = r"$\frac{\mathcal{H}'(x)}{\mathcal{H}(x)}$"
    ddH_label = r"$\frac{\mathcal{H}''(x)}{\mathcal{H}(x)}$"


    plot.compare_dH_and_ddH_over_H(x, dHp/Hp, ddHp/Hp, dH_label, ddH_label, m_dom, L_dom, 
                                   title=None, save=SAVE, temp=TEMP)


def eta_t_plot():
    ### Plot eta and t ### 
    eta, t = load_eta_and_t()
    eta_c = (eta / c).to(u.Gyr)

    mr_eq, mL_eq = equality_times()

    plot.plot_t_and_eta(x, t, eta_c, fname="t_and_eta_c.pdf", mr_eq=mr_eq, mL_eq=mL_eq, 
                        acc=None, save=SAVE, temp=TEMP)


def luminosity_distance(data):
    ### Return z and dL for a given data file
    #  limited to region with observational data 
    x_sim =  data[0]
    dL_sim   = (data[-1]*u.m).to(u.Gpc)

    z_sim = z_of_x(x_sim) 
    z_sim_mask = (z_sim <= 1.31) & (z_sim >= 0.005)
    z_sim = z_sim[z_sim_mask]

    dL_sim = dL_sim[z_sim_mask]
    
    return z_sim, dL_sim 


def plot_dL(plot_fit=False):
    ### Plot dL from data and compare with simulation with Planck parameters  
    ### plot_fit=True: Plot dL from the simulation with parameters obtained from supernovafit.
    #       Planck results included for comparison
    z, dL, dL_error = plot.load("supernovadata.txt", skiprows=1)
    planck_data     = plot.load("cosmology_dL.txt" , skiprows=1)
    fit_data        = plot.load("bestfit_cosmology_dL.txt" , skiprows=1)

    z_planck, dL_planck = luminosity_distance(planck_data)
    z_fit, dL_fit = luminosity_distance(fit_data)

    data    = [z,        dL, dL_error]
    planck  = [z_planck, dL_planck.value]
    fit     = [z_fit,    dL_fit.value]


    if plot_fit:
        plot.plot_dL(data, planck, fit, fname="dL_z_compare_fitted.pdf", save=SAVE, temp=TEMP)
    else:
        plot.plot_dL(data, planck, fname='dL_z_compare_planck.pdf', save=SAVE, temp=TEMP)




def load_supernovafit(burn):
    # Load result from supernova fit 
    # Leave out the first N=burn results  
    supernova_mcmc_results = plot.load("supernovafit.txt", skiprows=1+burn)
    chi2, h, OmegaM, OmegaK = supernova_mcmc_results
    return chi2, h, OmegaM, OmegaK


def supernova_fit_omegas(burn=1000):
    # Plot OmegaM-OmegaK confidence region.
    chi2, h, OmegaM, OmegaK = load_supernovafit(burn)
   
    OmegaLambda = 1 - OmegaM - OmegaK 
    chi2min = np.min(chi2)
    chi2_1sigma = chi2 < chi2min + 3.53
    chi2_2sigma = chi2 < chi2min + 8.02

    plot.plot_OmegaM_OmegaLambda_plane(OmegaM, OmegaLambda, 
                                       chi2_1sigma, chi2_2sigma, chi2_min=np.argmin(chi2),
                                       fname=f"mcmc_supernova_fit_Nburn{burn}.pdf", 
                                       save=SAVE, temp=TEMP)


def supernova_fit_H0_pdf(burn=1000):
    # Plot H0 pdf
    chi2, h, OmegaM, OmegaK = load_supernovafit(burn)

    H0 = h * (100*u.km / u.s / u.Mpc)

    H0_mean = np.mean(H0)
    H0_std  = np.std(H0) 
    H0_var  = np.var(H0) 

    print(f"H0 mean: {H0_mean:.5f}")
    print(f"H0 std : {H0_std:.5f}")

    
    H0_min = H0_mean - 4*H0_std 
    H0_max = H0_mean + 4*H0_std 
    N_bins = 100 
    bins   = np.linspace(H0_min, H0_max, N_bins)

    H0_gaussian_distr = (2*np.pi*H0_std**2)**(-1/2) * np.exp(-(bins - H0_mean)**2 / (2 * H0_var))

    plot.plot_H0_posterior_pdf(H0, bins, H0_gaussian_distr,
                                fname=f'H0_pdf_Nburn{burn}.pdf', save=SAVE, temp=TEMP)


def table():
    # make table. 
    # Used pandas initially, but changed the table so much that it's not applicable anymore.

    data("cosmology_times.txt") # load data with higher resolution
    MR_eq_idx, ML_eq_idx = equality_times(idx=True)
    acc_onset_idx = acceleration_onset(idx=True)
    x_today_idx = np.abs(x).argmin()

    eta, t = load_eta_and_t() 

    x_mr = x[MR_eq_idx]
    x_ml = x[ML_eq_idx]
    x_ac = x[acc_onset_idx]

    t_mr = t[MR_eq_idx]
    t_ml = t[ML_eq_idx]
    t_ac = t[acc_onset_idx]

    z_mr = z_of_x(x_mr)
    z_ml = z_of_x(x_ml)
    z_ac = z_of_x(x_ac)

    t_today = t[x_today_idx]
    eta_over_c_today = (eta[x_today_idx] / c).to(u.Gyr) 

    print(f"mr: x={x_mr:.3f}, z={z_mr:.3f}, t={t_mr.to(u.Gyr)/10**(-5):.3f}")
    print(f"ac: x={x_ac:.3f}, z={z_ac:.3f}, t={t_ac:.3f}")
    print(f"ml: x={x_ml:.3f}, z={z_ml:.3f}, t={t_ml:.3f}")
    print(f"t today= {t_today:.3f}")
    print(f"eta0   = {eta_over_c_today:.3f}")

    mr_eq = [x_mr, z_mr, t_mr.to(u.Gyr).value]
    ml_eq = [x_ml, z_ml, t_ml.value]
    acc   = [x_ac, z_ac, t_ac.value]

    #### Previous table generating with pandas     
    # print('making table')
    # plot.time_table(mr_eq, ml_eq, acc, t_today, eta_over_c_today, show=True, save=False)

    # Restore data set to the original one 
    data()


SAVE=True 
TEMP=True 
dH_ddH_over_H()
Hp_plot()
eta_plot()
eta_H_plot()
eta_t_plot()

plot_omegas()

plot_dL()
plot_dL(plot_fit=True)

supernova_fit_omegas()
supernova_fit_H0_pdf()

# table()
