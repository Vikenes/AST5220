import numpy as np 
import matplotlib.pyplot as plt 
import astropy.units as u 
from astropy.constants import c  
import os 
import plot

import warnings 
warnings.filterwarnings("ignore", category=DeprecationWarning )

save = False

simulation_data = plot.load("cosmology.txt")
x = simulation_data[0]



def load_H_parameters(convert_units=True):

    H_units = 100*u.km / u.s / u.Mpc 
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


def domination_eras():
    OmegaM, OmegaRel, OmegaLambda = load_omegas()

    matter_dom_onset_idx = np.where(OmegaM > OmegaRel+OmegaLambda)[0][0]
    Lambda_dom_onset_idx = np.where(OmegaLambda > OmegaM+OmegaRel)[0][0]

    matter_dom_onset = x[matter_dom_onset_idx]
    Lambda_dom_onset = x[Lambda_dom_onset_idx]

    return matter_dom_onset, Lambda_dom_onset



def Hp_plot(save):
    Hp = load_H_parameters()[0]
    title = r'$\mathcal{H}(x)\: \Big( \frac{100\,\mathrm{km/s}}{\mathrm{Mpc}} \Big) $'
    plot.plot_single_param(x, Hp, "compare_Hp.pdf", xlabel=r"x", 
                            xlim=[-12,0], ylim=[1e-1, 1e3], 
                            title=title, save=save)



def eta_plot(save):
    eta = load_eta_and_t()[0]
    title = r'$\eta(x)\: ( \mathrm{Mpc} ) $'

    plot.plot_single_param(x, eta, "compare_eta.pdf", xlabel=r"x", 
                            xlim=[-12,0], ylim=[1e0, 2e4], 
                            title=title, save=save)


def eta_H_plot(save):
    Hp = load_H_parameters()[0]
    eta = load_eta_and_t()[0]
    title = r'$\eta(x) \mathcal{H}(x) / c $'

    eta_Hp_over_c = (eta*Hp/c).to(1)
    mr_eq, mL_eq = domination_eras()
    plot.plot_single_param(x, eta_Hp_over_c, "compare_eta_H_over_c.pdf", 
                            mr_eq=mr_eq, mL_eq=mL_eq, legend=True,
                            xlabel=r"x", xlim=[-20,1], ylim=[0.75, 3], 
                            log=False, title=title, save=save)


def plot_omegas(save):
    OmegaM, OmegaRel, OmegaLambda = load_omegas()
    title = r"$\Omega_i(x)$"

    plot.plot_omega_params(x, OmegaM, OmegaRel, OmegaLambda, xlabel=r'$x$', 
                        fname="omega_i_of_x.pdf", title=title, save=save)
    


def dH_ddH_over_H(save):
    Hp, dHp, ddHp = load_H_parameters()
    dH_label = r"$\frac{1}{\mathcal{H}(x)} \frac{\mathrm{d} \mathcal{H}(x)}{\mathrm{d}x}$"
    m_dom, L_dom = domination_eras()

    plot.compare_dH_over_H(x, dHp/Hp, dH_label, m_dom, L_dom, title='tbd', save=save)

    ddH_label = r"$\frac{1}{\mathcal{H}(x)} \frac{\mathrm{d}^2 \mathcal{H}(x)}{\mathrm{d}x^2}$"
    plot.compare_ddH_over_H(x, ddHp/Hp, ddH_label, m_dom, L_dom, title='tbd', save=save)


def eta_t_plot(save):
    eta, t = load_eta_and_t()
    eta_c = (eta / c).to(u.Gyr)


    plot.plot_t_and_eta(x, t, eta_c, fname="t_and_eta_c.pdf", save=save)



def luminosity_distance(save):
    z, dL, dL_error = plot.load("supernovadata.txt", skiprows=1)
    supernova_x     = plot.load("cosmology_dL.txt" , skiprows=1)
    x_sim_dL = supernova_x[0]
    dL_sim = (supernova_x[-1]*u.m).to(u.Gpc)

    z_sim = np.exp(-x_sim_dL) - 1 
    z_sim_mask = (z_sim <= 1.31) & (z_sim >= 0.005)
    z_sim = z_sim[z_sim_mask]
    dL_sim = dL_sim[z_sim_mask]
    

    plot.plot_dL(z, dL, dL_error, z_sim, dL_sim.value, fname="dL_z_compare_log.pdf", save=save)


def load_supernovafit(burn):
    # supernova_mcmc_results = plot.load("supernovafit_h0_new.txt", skiprows=1)
    supernova_mcmc_results = plot.load("supernovafit.txt", skiprows=1)
    chi2, h, OmegaM, OmegaK = supernova_mcmc_results[:,burn:]
    return chi2, h, OmegaM, OmegaK


def supernova_fit_omegas(save, burn=1000):

    chi2, h, OmegaM, OmegaK = load_supernovafit(burn)
    OmegaLambda = 1 - OmegaM - OmegaK 
    chi2min = np.min(chi2)
    chi2_1sigma = chi2 < chi2min + 3.53
    chi2_2sigma = chi2 < chi2min + 8.02

    plot.plot_OmegaM_OmegaLambda_plane(OmegaM, OmegaLambda, chi2_1sigma, chi2_2sigma, 
                    fname=f"mcmc_supernova_fit_Nburn{burn}.pdf", save=save)


def supernova_fit_H0_pdf(save, burn=1000):
    chi2, h, OmegaM, OmegaK = load_supernovafit(burn)

    H0 = h * (100*u.km / u.s / u.Mpc)

    H0_mean = np.mean(H0)
    H0_std  = np.std(H0) 
    H0_var  = np.var(H0) 
    
    H0_min = H0_mean - 4*H0_std 
    H0_max = H0_mean + 4*H0_std 
    N_bins = 100 
    bins   = np.linspace(H0_min, H0_max, N_bins)

    H0_gaussian_distr = (2*np.pi*H0_std**2)**(-1/2) * np.exp(-(bins - H0_mean)**2 / (2 * H0_var))

    plot.plot_H0_posterior_pdf(H0, bins, H0_gaussian_distr,
                                fname=f'H0_pdf_Nburn{burn}.pdf', save=save)


    # s2 = np.var(OmegaLambda)
    # mu = np.mean(OmegaLambda)

    # plt.hist(OmegaLambda[burn:], bins=40, density=True)
    # x = np.linspace(0, 1.5, 100)
    # plt.plot(x, 1/np.sqrt(2*np.pi*s2) * np.exp(- (x-mu)**2 / (2 * s2)))
    # plt.vlines(0.685, 0, 10, color='red')
    # plt.ylim(0, 2.6)
    # plt.xlim(-0.1, 1.6)
    # plt.show()





# luminosity_distance(save)
# exit()
# dH_ddH_over_H(save)
# Hp_plot(save)
# eta_plot(save)
# eta_H_plot(save)
# eta_t_plot(save)
# plot_omegas(save)


# supernova_fit_omegas(save)
# supernova_fit_H0_pdf(save)

