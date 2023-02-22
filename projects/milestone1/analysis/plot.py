import numpy as np 
import matplotlib.pyplot as plt 
import astropy.units as u 
from astropy.constants import c  
import os 

plt.style.use("seaborn")

# Paths 
here = os.path.abspath(".")
# data_path = here + "/../output/data/"
# latex_path = here + "/../../latex/"
# temp_path = here + "/../../output/plots/temp/"
# plot_path = here +"/../../output/plots/pdfs/"


#   rc and plot params
TICKLABELSIZE = 25
LABELSIZE = 25
LEGENDSIZE = 20
TITLESIZE = 30

#   Set rc-params
plt.rc("legend", fontsize=LEGENDSIZE, fancybox=True, loc="best", frameon=True, edgecolor="black")
plt.rc("font", size=25)
plt.rc("axes", titlesize=TITLESIZE, labelsize=LABELSIZE)
plt.rc("xtick", labelsize=TICKLABELSIZE)
plt.rc("ytick", labelsize=TICKLABELSIZE)
# plt.rc("tickparams")
plt.rc("lines", linewidth=2.5)
plt.rc('text', usetex=True)
# plt.rc("label", fontsize=20)

plt.rcParams['savefig.bbox'] = 'tight'
plt.rcParams['font.family'] = 'Times New Roman'


def save_push(fig, pdf_name, save=True, push=False, show=False, tight=True):
    """
    This function handles whether you want to show,
    save and/or push the file to git.
    Args:
        fig (matplotlib.figure): Figure you want to handle
        pdfname (string): Name of output pdf file
        args (argparse)
    """
    if tight:
        fig.tight_layout()
    pdfname = pdf_name.replace('.pdf', '').strip() + ".pdf"
    file = plot_path + pdfname
    if save:
        print(f'Saving plot: {file}')
        fig.savefig(file)
        if png_duplicate:
            png_fname = temp_path + pdfname.replace('.pdf', '.png')
            fig.savefig(png_fname)
            os.system(f"git add {png_fname}")
    if push:
        os.system(f"git add {file}")
        os.system("git commit -m 'upload plot'")
        os.system("git push")
    if show:
        plt.show()
    if not show:
        plt.close()
    else:
        plt.close()



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