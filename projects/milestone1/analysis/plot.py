import numpy as np 
import matplotlib.pyplot as plt 
import astropy.units as u 
from astropy.constants import c  
import os 
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter


import warnings 
warnings.filterwarnings("ignore", category=DeprecationWarning )

### Style and formatting of plots 
plt.style.use("seaborn")



#   rc and plot params
TICKLABELSIZE = 20
LABELSIZE     = 25
LEGENDSIZE    = 20
TITLESIZE     = 30
LABELPAD      = 12

#   Set rc-params
plt.rc("legend", fontsize=LEGENDSIZE, fancybox=True, loc="best", frameon=True, edgecolor="black")
plt.rc("font", size=25)
plt.rc("axes", titlesize=TITLESIZE, labelsize=LABELSIZE, labelpad=LABELPAD)
plt.rc("xtick", labelsize=TICKLABELSIZE)
plt.rc("ytick", labelsize=TICKLABELSIZE)
plt.rc("ytick.minor", pad=10)
plt.rc("lines", linewidth=2.5)
plt.rc('text', usetex=True)
# plt.rc("label", fontsize=20)

plt.rcParams['savefig.bbox'] = 'tight'
plt.rcParams['font.family'] = 'Times New Roman'



# Paths 
here = os.path.abspath(".")
data_path = here + "/../data/"
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
        print(f'Saving plot: {pdfname}')
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



def plot_single_param(x, quantity, fname, mr_eq=None, mL_eq=None,
                        xlabel=None, ylabel=None, title=None, legend=False,
                        xlim=None, ylim=None, log=True, save=True, push=False):

    fig, ax = plt.subplots(figsize=(10,8))
    ax.plot(x, quantity)

    if ylim is None:
        ymin = np.min(quantity)
        ymax = np.max(quantity)
    else:
        ymin, ymax = ylim


    if mr_eq is not None:
        ax.vlines(mr_eq, ymin, ymax, ls='dashed',
                    color='red', label=r'$\Omega_\mathrm{rel}=\Omega_m$')

    if mL_eq is not None:
        ax.vlines(mL_eq, ymin, ymax, ls='dashed',
                    color='green', label=r'$\Omega_m=\Omega_\Lambda$')

    ax.set_ylim(ylim)
    ax.set_xlim(xlim)

    if log:
        ax.set_yscale('log')
    
    set_ax_info(ax, xlabel, ylabel, title, legend=legend)
    
    save_push(fig, fname, save, push)
    

def plot_omega_params(x, m, r, L, fname, xlabel=None, ylabel=None, title=None, 
                    xlim=[-20,5], ylim=[0,1.1], save=True, push=False):

    fig, ax = plt.subplots(figsize=(10,8))

    ax.plot(x, m, label=r'$\Omega_\mathrm{m} = \Omega_\mathrm{b} + \Omega_\mathrm{CDM}$')
    ax.plot(x, r, label=r'$\Omega_\mathrm{rel} = \Omega_\gamma + \Omega_\nu$')
    ax.plot(x, L, label=r'$\Omega_\Lambda$')



    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    set_ax_info(ax, xlabel, ylabel, title, legend=True)
    
    save_push(fig, fname, save, push)



def compare_dH_over_H(x, dH_H, H_label, x_mr_eq, x_mL_eq, 
                        title=None, xlim=[-20,5], save=True, push=False):
    
    fig, ax = plt.subplots(figsize=(12,10))

    ax.plot(x, dH_H, label=H_label)
    ax.hlines(-1, xmin=xlim[0], xmax=x_mr_eq,  ls='dashed',color='r', label=r'$\omega=1/3$')
    ax.hlines(-1/2, xmin=x_mr_eq, xmax=x_mL_eq,ls='dashed',color='green', label=r'$\omega=0$')
    ax.hlines(1, xmin=x_mL_eq, xmax=xlim[-1],  ls='dashed',color='orange', label=r'$\omega=-1$')

    ax.set_xlim(xlim)

    set_ax_info(ax, xlabel='$x$', title=title)
    save_push(fig, 'dH_over_H.pdf', save)


        
def compare_ddH_over_H(x, ddH_H, H_label, x_mr_eq, x_mL_eq, 
                        title=None, xlim=[-20,5], save=True, push=False):
    
    fig, ax = plt.subplots(figsize=(12,10))

    ax.plot(x, ddH_H, label=H_label)
    ax.hlines(1, xmin=xlim[0], xmax=x_mr_eq,  ls='dashed',color='r', label=r'$\omega=1/3$')
    ax.hlines(1/4, xmin=x_mr_eq, xmax=x_mL_eq,ls='dashed',color='green', label=r'$\omega=0$')
    ax.hlines(1, xmin=x_mL_eq, xmax=xlim[-1],  ls='dashed',color='orange', label=r'$\omega=-1$')

    ax.set_xlim(xlim)

    set_ax_info(ax, xlabel='$x$', title=title)
    save_push(fig, 'ddH_over_H.pdf', save)


def plot_t_and_eta(x, t, etac, fname, xlim=[-20,5], save=True):
    fig, ax = plt.subplots(figsize=(8,6))

    ax.set_xlim(xlim)
    ax.plot(x, t, label=r'$t(x)$')
    ax.plot(x, etac, label=r'$\eta(x)/c$')
    
    ylabel=r'Time $[\mathrm{Gyr}]$'
    xlabel=r"$x$"


    set_ax_info(ax, xlabel, ylabel)
    save_push(fig, fname, save=save)



def plot_dL(z_data, dL_data, dL_error, z_sim, dL_sim, 
            fname="dL_z_compare.pdf", save=True):
    fig, ax = plt.subplots(figsize=(10,8))


    ax.plot(z_sim, dL_sim/z_sim, label='Simulation')
    ax.errorbar(z_data, dL_data/z_data, yerr=dL_error/z_data, barsabove=True, fmt='x', 
                capthick=1.5, capsize=5, elinewidth=2, color='r',
                label='Data')

    ax.set_xscale('log')
    ax.set_yscale('log')
    xlabel=r'$z$'
    ylabel=r'$d_L(z)/z \:[\mathrm{Gpc}]$'

    # for axis in [ax.xaxis, ax.yaxis]:

    set_ax_info(ax, xlabel, ylabel, title="Supernova fit")
    # formatter = ScalarFormatter()
    # formatter.set_scientific(False)
    ax.yaxis.set_minor_formatter(ScalarFormatter())
    # ax.yaxis.set_minor_formatter(FormatStrFormatter("%d"))
    # ax.yaxis.set_major_formatter(ScalarFormatter())

    # plt.ticklabel_format(axis='y', style='plain')
    save_push(fig, fname, save=save)



def plot_OmegaM_OmegaLambda_plane(OmegaM, OmegaLambda, chi2_1sigma, chi2_2sigma, fname, save=True):

    fig, ax = plt.subplots(figsize=(12,10))

    ax.plot(OmegaM[chi2_2sigma], OmegaLambda[chi2_2sigma],'ro',ms=2)
    ax.plot(OmegaM[chi2_1sigma], OmegaLambda[chi2_1sigma],'bo',ms=2)
    ax.plot([], 'ro', label=r'$2\sigma$')
    ax.plot([], 'bo', label=r'$1\sigma$')
    omega_k_zero = np.linspace(0, 1, 20)
    ax.plot(omega_k_zero, 1 - omega_k_zero, 'k--', label=r'$\Omega_k=0$')

    ax.set_xlim(0,1)
    ax.set_ylim(0,1.5)
    set_ax_info(ax, xlabel=r"$\Omega_m$", ylabel=r"$\Omega_\Lambda$")
    save_push(fig, fname, save)


def plot_H0_posterior_pdf(H0, bins, H0_gaussian, fname, save=True):
    fig, ax = plt.subplots(figsize=(12,8))

    ax.plot(bins, H0_gaussian, label=r"$H_0\sim \mathcal{N}(\mu,\sigma^2)$")
    ax.hist(H0, bins=bins, density=True, edgecolor='k')
    set_ax_info(ax, xlabel=r"$H_0\:[\mathrm{km/s/Mpc}]$")
    save_push(fig, fname, save)



if __name__=='__main__':
    pass 