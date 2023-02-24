import numpy as np 
import matplotlib.pyplot as plt 
import astropy.units as u 
from astropy.constants import c  
import os 
import plot

import warnings 
warnings.filterwarnings("ignore", category=DeprecationWarning )


simulation_data = plot.load("cosmology_large.txt")
x = simulation_data[0]
print(len(x), x[0], x[-1])
# exit()
# def x(file="cosmology15.txt"):
    # return plot.load(file)[0]

def load_H_parameters(convert_units=True):
    # simulation_results = plot.load(file)
    # x = simulation_results[0]

    H_units = 100*u.km / u.s / u.Mpc 
    H_params = simulation_data[2:5] / u.s 

    if convert_units: 
        H_params = H_params.to(100*u.km / u.s / u.Mpc) 

    Hp, dHp_dx, ddHp_ddx = H_params 

    return Hp, dHp_dx, ddHp_ddx


def load_eta_and_t(convert_units=True):
    eta = simulation_data[1] * u.m 
    t   = simulation_data[-1] * u.s 

    if convert_units:
        eta = eta.to(u.Mpc)
        t = t.to(u.Gyr)

    return eta, t 


def load_omegas():
    OmegaB, OmegaCDM, OmegaLambda = simulation_data[5:8]
    OmegaR, OmegaNu, OmegaK = simulation_data[8:11]

    OmegaM = OmegaB + OmegaCDM
    OmegaRel = OmegaR + OmegaNu 
    return OmegaM, OmegaRel, OmegaLambda 


def Hp_plot(save=True):
    Hp = load_H_parameters()[0]
    title = r'$\mathcal{H}(x)\: \Big( \frac{100\,\mathrm{km/s}}{\mathrm{Mpc}} \Big) $'

    plot.plot_single_param(x, Hp, "compare_Hp.pdf", xlabel=r"x", 
                            xlim=[-12,0], ylim=[1e-1, 1e3], 
                            title=title, save=save)



def eta_plot(save=True):
    eta = load_eta_and_t()[0]
    title = r'$\eta(x)\: ( \mathrm{Mpc} ) $'

    plot.plot_single_param(x, eta, "compare_eta.pdf", xlabel=r"x", 
                            xlim=[-12,0], ylim=[1e0, 2e4], 
                            title=title, save=save)


def eta_H_plot(save=True):
    Hp = load_H_parameters()[0]
    eta = load_eta_and_t()[0]
    title = r'$\eta(x) \mathcal{H}(x) / c $'

    eta_Hp_over_c = (eta*Hp/c).to(1)

    plot.plot_single_param(x, eta_Hp_over_c, "compare_eta_H_over_c.pdf", xlabel=r"x", 
                            xlim=[-20,0], ylim=[0.75, 3], log=False,
                            title=title, save=save)


def plot_omegas(save=True):
    OmegaM, OmegaRel, OmegaLambda = load_omegas()
    title = r"$\Omega_i(x)$"

    plot.plot_omega_params(x, OmegaM, OmegaRel, OmegaLambda, xlabel=r'$x$', 
                        fname="omega_i_of_x.pdf", title=title, save=save)
    

# Hp_plot()
# eta_plot()
eta_H_plot(False)
# plot_omegas()