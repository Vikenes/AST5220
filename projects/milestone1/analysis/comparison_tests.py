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
    # plt.plot(x, dHp/Hp)
    # plt.vlines(m_dom, -1.2, 1.2)
    # plt.vlines(L_dom, -1.2, 1.2)

    # plt.show()
    ddH_label = r"$\frac{1}{\mathcal{H}(x)} \frac{\mathrm{d}^2 \mathcal{H}(x)}{\mathrm{d}x^2}$"

    plot.compare_ddH_over_H(x, ddHp/Hp, ddH_label, m_dom, L_dom, title='tbd', save=save)


def eta_t_plot(save):
    eta, t = load_eta_and_t()
    eta_c = (eta / c).to(u.Gyr)

    plot.plot_t_and_eta(x, t, eta_c, fname="t_and_eta_c.pdf", save=save)



def luminosity_distance(save):
    z, dL, dL_error = np.loadtxt("../data/supernovadata.txt", unpack=True, skiprows=1)
    x_min = np.min(-np.log(1+z))
    x_max = np.max(-np.log(1+z))
    x_mask = (x >= x_min) & (x <= x_max)
    z_sim = np.exp(-x[x_mask]) - 1

    x0_idx = np.argmin(np.abs(x))
    eta = load_eta_and_t()[0].to(u.Gpc)
    eta0 = eta[x0_idx]
    etax = eta[x_mask]
    r = eta0 - etax 
    dL_sim = r / np.exp(x[x_mask])



    plot.plot_dL(z, dL, dL_error, z_sim, dL_sim, save)


luminosity_distance(True)

# dH_ddH_over_H(save)
# Hp_plot(save)
# eta_plot(save)
# eta_H_plot(save)
# eta_t_plot(save)
# plot_omegas(save)

# y_min = None 

# print(np.min([int(y_min), 1]))