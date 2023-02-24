import numpy as np 
import matplotlib.pyplot as plt 
import astropy.units as u 
from astropy.constants import c  
import os 


import warnings 
warnings.filterwarnings("ignore", category=DeprecationWarning )

### Style and formatting of plots 
plt.style.use("seaborn")



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



# Paths 
here = os.path.abspath(".")
data_path = here + "/../out_data/"
# latex_path = here + "/../../latex/"
# temp_path = here + "/../../output/plots/temp/"
fig_path = here +"/figures/"

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
    file = fig_path + pdfname
    
    if save:
        print(f'Saving plot: {file}')
        fig.savefig(file)
    else:
        plt.show()
    if push:
        os.system(f"git add {file}")
        os.system("git commit -m 'upload plot'")
        os.system("git push")

    plt.close()


def set_ax_info(ax, xlabel, ylabel=False, title=None, legend=True):
    """Write title and labels on an axis with the correct fontsizes.
    Args:
        ax (matplotlib.axis): the axis on which to display information
        title (str): the desired title on the axis
        xlabel (str): the desired lab on the x-axis
        ylabel (str): the desired lab on the y-axis
    """
    ax.set_xlabel(xlabel)
    if ylabel != False:
        ax.set_ylabel(ylabel)
    ax.set_title(title)
    # ax.tick_params(axis='both', which='major', labelsize=15)
    # ax.yaxis.get_offset_text().set_fontsize(15)
    # try:
        # ax.ticklabel_format(style=style)
    # except AttributeError:
        # pass
    if legend:
        ax.legend()




# -----------------------------------------------------------------------------
#   General plot code above
#   XXX
#   Specific plot code below
# -----------------------------------------------------------------------------



def load(file, folder=data_path, skiprows=0):
    # Tired of repeating unpack, delimiter, skiprows for all... 
    return np.loadtxt(folder + file, unpack=True, skiprows=skiprows)



def plot_single_param(x, quantity, fname, xlabel, ylabel=None, title=None, 
                        xlim=None, ylim=None, log=True, save=True, push=False):

    fig, ax = plt.subplots(figsize=(10,8))
    ax.plot(x, quantity)
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    if log:
        ax.set_yscale('log')
    
    set_ax_info(ax, xlabel, ylabel, title, legend=False)
    
    save_push(fig, fname, save, push)
    

def plot_omega_params(x, m, r, L, fname, xlabel=None, ylabel=None, title=None, 
                    xlim=[-20,5], ylim=[0,1.1], save=True, push=False):

    fig, ax = plt.subplots(figsize=(10,8))

    ax.plot(x, m, label=r'$\Omega_\mathrm{m} = \Omega_\mathrm{b} + \Omega_\mathrm{CDM}$')
    ax.plot(x, r, label=r'$\Omega_\mathrm{rel} = \Omega_\mathrm{r} + \Omega_\nu$')
    ax.plot(x, L, label=r'$\Omega_\Lambda$')



    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    set_ax_info(ax, xlabel, ylabel, title, legend=True)
    
    save_push(fig, fname, save, push)

def comparison_tests(param_test=[]):
    cosmology_data = load("cosmology.txt")
    x = cosmology_data[0]
    eta = (cosmology_data[1] * u.m).to(u.Mpc)
    Hp, dHp_dx = (cosmology_data[2:4] / u.s).to(100*u.km / u.s / u.Mpc)

    OmegaB, OmegaCDM, OmegaLambda, OmegaR = cosmology_data[4:8]
    OmegaNu, OmegaK = cosmology_data[8:-1]
    Omega_tot = OmegaNu + OmegaB + OmegaCDM + OmegaK + OmegaLambda + OmegaR

    t = (cosmology_data[-1]*u.s).to(u.Gyr)

    dhp_h = lambda w: - (1 + 3*w) / (2 * np.exp(x))

    plt.plot(x, dHp_dx/Hp)
    # plt.plot(x, dhp_h(0), '--', label='m')
    # plt.plot(x, dhp_h(-1), '--', label=r'$\Lambda$')
    # plt.plot(x, dhp_h(1/3), '--', label='rad')

    # plt.yscale('log')
    plt.legend()
    # plt.xlim(-5,0.2)
    # plt.ylim(-10,10)
    plt.show()







def supernova_fit():
    chi2, h, OmegaM, OmegaK = load("results_supernovafitting_1e4.txt", skiprows=1)
    OmegaLambda = 1 - OmegaM - OmegaK 
    chi2min = np.min(chi2)
    chi2_1sigma = chi2 < chi2min + 3.53
    chi2_2sigma = chi2 < chi2min + 8.02 
    # flat = OmegaM == OmegaLambda
    # flat = np.abs(OmegaM + OmegaLambda).argmin()
    plt.plot(OmegaM[chi2_2sigma], OmegaLambda[chi2_2sigma],'ro',ms=3, label=r'$2\sigma$')
    plt.plot(OmegaM[chi2_1sigma], OmegaLambda[chi2_1sigma],'bo',ms=3, label=r'$1\sigma$')
    plt.plot(OmegaM, 1 - OmegaM, 'k--', label='flat')
    # plt.plot(OmegaM, 1 - OmegaM, label='meh')

    plt.xlim(0,1)
    plt.ylim(0,1.5)
    plt.legend()
    plt.show()

# comparison_tests()
# supernova_fit()

"""

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
"""

if __name__=='__main__':
    pass 