import numpy as np 
import matplotlib.pyplot as plt 
import astropy.units as u 
from astropy.constants import c  
import os 

plt.style.use("seaborn")

cosmology_data = np.loadtxt("../cosmology.txt", unpack=True)
x = cosmology_data[0]
eta = (cosmology_data[1] * u.m).to(u.Mpc)
Hp, dHp_dx = (cosmology_data[2:4] / u.s).to(100*u.km / u.s / u.Mpc)

OmegaB, OmegaCDM, OmegaLambda, OmegaR = cosmology_data[4:8]
OmegaNu, OmegaK = cosmology_data[8:-1]
Omega_tot = OmegaNu + OmegaB + OmegaCDM + OmegaK + OmegaLambda + OmegaR

t = (cosmology_data[-1]*u.s).to(u.Gyr)

def plot_Hp():
    plt.plot(x, Hp)
    plt.xlim(-12,0)
    plt.ylim(1e-1, 1e3)
    plt.yscale('log')
    plt.show()

def plot_eta():
    plt.plot(x, eta)
    # plt.xlim(-12,0)
    # plt.ylim(1, 2e4)
    plt.yscale('log')
    plt.show()

def plot_etaH_c():
    # print(((eta*Hp/c).decompose())[0])
    plt.plot(x, (eta*Hp/c).to(1))
    plt.xlim(-12,0)
    # plt.ylim(1, 2e4)
    # plt.yscale('log')
    plt.show()

def plot_omegas():
    plt.plot(x, OmegaR+OmegaNu, label=r"$\Omega_r$")
    plt.plot(x, OmegaB+OmegaCDM, label=r"$\Omega_m=\Omega_b+\Omega_{CDM}$")
    plt.plot(x, OmegaLambda, label=r"$\Omega_\Lambda$")
    plt.plot(x, Omega_tot,'k--')
    plt.legend()
    plt.show()

def plot_t():
    plt.plot(x, t)
    plt.yscale('log')
    print(f't(0) = {t[-1]:.3f}')
    plt.show()

# plot_Hp()
# plot_eta()
# plot_etaH_c()
# plot_omegas()
plot_t()